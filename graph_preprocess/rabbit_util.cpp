//
// Created by atrostan on 07/09/22.
//

#include "rabbit_util.h"


namespace rabbit_order {
	namespace aux {

//==============================================================================
//
// UTILITIES
//
//==============================================================================

//------------------------------------------------------------------------------
//
// Helpers of `RO_DIE`
//
//------------------------------------------------------------------------------

		struct die_t {
			die_t(const char *file, uint64_t line, const char *func) {
				std::cerr << file << ':' << line << '(' << func << ") [FATAL ERROR] ";
			}

			~die_t() {
				std::cerr << std::endl;
				exit(EXIT_FAILURE);
			}
		};

		template<typename T>
		die_t operator<<(die_t d, const T &x) {
			std::cerr << x;
			return d;
		}
	}
}

off_t file_size(const int fd) {
	struct stat st;
	if (fstat(fd, &st) != 0)
		RO_DIE << "stat(2): " << strerror(errno);
	return st.st_size;
}

struct file_desc {
	int fd;

	file_desc(const std::string &path) {
		fd = open(path.c_str(), O_RDONLY);
		if (fd == -1)
			RO_DIE << "open(2): " << strerror(errno);
	}

	~file_desc() {
		if (close(fd) != 0)
			RO_DIE << "close(2): " << strerror(errno);
	}
};

struct mmapped_file {
	size_t size;
	const char *data;

	mmapped_file(const file_desc &fd) {
		size = file_size(fd.fd);
		data = static_cast<char *>(
				mmap(NULL, size, PROT_READ, MAP_PRIVATE, fd.fd, 0));
		if (data == NULL)
			RO_DIE << "mmap(2): " << strerror(errno);
	}

	~mmapped_file() {
		if (munmap(const_cast<char *>(data), size) != 0)
			RO_DIE << "munmap(2): " << strerror(errno);
	}
};

struct edge_parser {
	const char *const strfirst;
	const char *const strlast;
	const char *crr;

	edge_parser(const char *const first, const char *const last)
			: strfirst(first), strlast(last), crr(first) {}

	template<typename OutputIt>
	void operator()(OutputIt dst) {
		while (crr < strlast) {
			eat_empty_lines();
			if (crr < strlast)
				*dst++ = eat_edge();
		}
	}

	std::pair<ul, ul> eat_edge() {
		const ul s = eat_id();
		eat_separator();
		const ul t = eat_id();
		return std::make_pair(s, t);  // FIXME: edge weight is not supported so far
	}

	ul eat_id() {
		//
		// Naive implementation is faster than library functions such as `atoi` and
		// `strtol`
		//
		const auto _crr = crr;
		ul v = 0;
		for (; crr < strlast && std::isdigit(*crr); ++crr) {
			const ul _v = v * 10 + (*crr - '0');
			if (_v < v)  // overflowed
				RO_DIE << "Too large vertex ID at line " << crr_line();
			v = _v;
		}
		if (_crr == crr) {  // If any character has not been eaten
			RO_DIE << "Invalid vertex ID at line " << crr_line();
		}
		return v;
	}

	void eat_empty_lines() {
		while (crr < strlast) {
			if (*crr == '\r') ++crr;                                // Empty line
			else if (*crr == '\n') ++crr;                                // Empty line
			else if (*crr == '#') crr = std::find(crr, strlast, '\n');  // Comment
			else if (*crr == '%') crr = std::find(crr, strlast, '\n');  // Comment
			else break;
		}
	}

	void eat_separator() {
		while (crr < strlast && (*crr == '\t' || *crr == ',' || *crr == ' '))
			++crr;
	}

	// Only for error messages
	size_t crr_line() {
		return std::count(strfirst, crr, '\n');
	}
};

/*
 * par_read thanks to rabbit order repo https://github.com/araij/rabbit_order
 */
std::vector<std::pair<ul, ul>> par_read(const std::string &path) {
	const file_desc fd(path);
	const mmapped_file mm(fd);
	const int nthread = omp_get_max_threads();
	const size_t zchunk = 1024 * 1024 * 64;  // 64MiB
	const size_t nchunk = mm.size / zchunk + (mm.size % zchunk > 0);

	//
	// For load balancing, partition the file into small chunks (whose size is
	// defined as `zchunk`) and dynamically assign the chunks into threads
	//

	std::vector<std::deque<std::pair<ul, ul>>> eparts(nthread);
#pragma omp parallel for schedule(dynamic, 1)
	for (size_t i = 0; i < nchunk; ++i) {
		const char *p = mm.data + zchunk * i;
		const char *q = mm.data + std::min(zchunk * (i + 1), mm.size);

		// Advance pointer `p` to the end of a line because it is possibly at the
		// middle of the line
		if (i > 0) p = std::find(p, q, '\n');

		if (p < q) {  // If `p == q`, do nothing
			q = std::find(q, mm.data + mm.size, '\n');  // Advance `q` likewise
			edge_parser(p, q)(std::back_inserter(eparts[omp_get_thread_num()]));
		}
	}

	// Compute indices to copy each element of `eparts` to
	std::vector<size_t> eheads(nthread + 1);
	for (int t = 0; t < nthread; ++t)
		eheads[t + 1] = eheads[t] + eparts[t].size();

	// Gather the edges read by each thread to a single array
	std::vector<std::pair<ul, ul>> edges(eheads.back());
#pragma omp parallel for schedule(guided, 1)
	for (int t = 0; t < nthread; ++t)
		boost::copy(eparts[t], edges.begin() + eheads[t]);

	return edges;
}