//
// Created by atrostan on 14/09/22.
//


/**
* given a graph that has been reordered using slashburn, compute
 * the fur-hilbert order of the edges contained within it's k*i core
 * (upper, left square in the adjacency matrix)
*/

#include "sb_furhilbert.h"


int main(int argc, char *argv[]) {
	opterr = 0;
	int opt;
	uint32_t n = 0;
	while ((opt = getopt(argc, argv, "n:")) != -1) {
		switch (opt) {
			case 'n':
				n = atoi(optarg);
				break;
		}
	}
	uint32_t imin;
	uint32_t imax;
	uint32_t jmin;
	uint32_t jmax;
	uint32_t i, j;

	uint64_t hilbert_counter = 0;
	i = j = 0;

//	int n_threads = omp_get_max_threads();
	int n_threads = 1;

	int n_quads = 1;
	int critical_depth = 0;
	while (n_quads < n_threads) {
		critical_depth += 1;
		n_quads = n_quads * 4;
	}
	fmt::print("n: {}\n", n);
	fmt::print("n_threads: {}\n", n_threads);

	std::vector<Quadrant> qs(n_quads);
	fmt::print("n_quads: {}\n", n_quads);

	QuadrantCalculator qc = QuadrantCalculator(critical_depth, n, qs);
	std::string debug_path = "/home/atrostan/Workspace/repos/congenial-enigma/graph_preprocess/cmake-build-debug/debug/hilbert";
	std::ofstream outfile(debug_path);

	for (auto &q: qs) {
		int i = 0;
		int j = 0;
		imin = q.start_x;
		jmin = q.start_y;
		imax = q.end_x;
		jmax = q.end_y;
		print_quad(q);
		FUR_HILBERT_FOR(i, j, imin, imax, jmin, jmax)
						{
							outfile << i << " " << j << "\n";
						}
		FUR_HILBERT_END(i, j);

	}
	outfile.close();
}