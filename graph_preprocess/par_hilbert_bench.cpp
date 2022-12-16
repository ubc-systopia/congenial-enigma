//
// Created by atrostan on 13/12/22.
//

#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include <getopt.h>
#include <cstdint>
#include <vector>
#include <string>
#include <omp.h>
#include <likwid-marker.h>
#include <likwid.h>
#include "fmt/core.h"
#include "fmt/ranges.h"
#include "io.h"
#include "omp.h"
#include "util.h"
#include <bitset>
#include <limits.h>


struct Quad {
	uint32_t qx;
	uint32_t qy;
	uint32_t q_idx;
	uint32_t nnz;
	std::vector<uint32_t> edges;
};

void write_binary_qs(std::string filename, std::vector<Quad> &qs, uint32_t q_side_len, uint32_t n_quads_per_side) {
	std::ofstream out(filename, std::ios::binary | std::ios::out | std::ios::trunc);
	// write the number of quadrants, sidelength of quadrants, and number of quads per side
	uint32_t n_quads = qs.size();
	out.write(reinterpret_cast<char *>(&n_quads), sizeof(uint32_t));
	out.write(reinterpret_cast<char *>(&q_side_len), sizeof(uint32_t));
	out.write(reinterpret_cast<char *>(&n_quads_per_side), sizeof(uint32_t));

	// write each quadrant and its metadata
	for (uint32_t i = 0; i < n_quads; ++i) {
		uint32_t qx = qs[i].qx;
		uint32_t qy = qs[i].qy;
		uint32_t q_idx = qs[i].q_idx;
		uint32_t nnz = qs[i].nnz;

		out.write(reinterpret_cast<char *>(&qx), sizeof(uint32_t));
		out.write(reinterpret_cast<char *>(&qy), sizeof(uint32_t));
		out.write(reinterpret_cast<char *>(&q_idx), sizeof(uint32_t));
		out.write(reinterpret_cast<char *>(&nnz), sizeof(uint32_t));

		// write the flattened edges in this quadrant
		out.write(
			reinterpret_cast<const char *>(qs[i].edges.data()),
			qs[i].nnz * sizeof(uint32_t) * 2);
	}
	out.close();
}

uint64_t highest_power_of_2_greater_than(uint32_t x) {
	// check for the set bits
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	// Then we remove all but the top bit by xor'ing the
	// string of 1's with that string of 1's shifted one to
	// the left, and we end up with just the one top bit
	// followed by 0's.
	x ^= (x >> 1);
	return x * 2;
}


ostream &operator<<(ostream &os, const Quad q) {
//	os << g.id << endl << g.scores.size() << endl;
	os << "Q" << q.q_idx << ": (" << q.qx << ", " << q.qy << "), "
	   << q.nnz;
	return os;
}

/**
 * For each quadrant in a vector of quadrant:
 * Reorder the edges contained in that quadrant using the hilbert order
 * @param qs
 * @param side_len
 */
void horder_edges_in_qs(std::vector<Quad> &qs, uint32_t side_len) {

#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < qs.size(); ++i) {
		// create a temporary vector to store the reordered edges
		uint32_t nnz = qs[i].nnz;
		struct IndexedEdge {
			uint32_t src;
			uint32_t dest;
			uint32_t h_idx;
		};
		std::vector<IndexedEdge> copy(nnz);
		for (uint32_t j = 0; j < nnz; ++j) {
			uint32_t src = qs[i].edges[j * 2];
			uint32_t dest = qs[i].edges[j * 2 + 1];
			copy[j].src = src;
			copy[j].dest = dest;
			copy[j].h_idx = xy2d(side_len, src, dest);
		}
		std::sort(
			copy.begin(),
			copy.end(),
			[](const IndexedEdge &a, const IndexedEdge &b) -> bool {
				return a.h_idx < b.h_idx;
			}
		);
		// copy the (now sorted) edges back to the edges array
		for (uint32_t j = 0; j < nnz; ++j) {
			qs[i].edges[j * 2] = copy[j].src;
			qs[i].edges[j * 2 + 1] = copy[j].dest;
		}
	}
}

void horder_qs(std::vector<Quad> &qs, uint32_t side_len) {
	uint64_t n = hyperceiling(side_len);
//	fmt::print("side_len: {}\n", side_len);
//	fmt::print("n: {}\n", n);
#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < qs.size(); ++i) {
		uint32_t qx = qs[i].qx;
		uint32_t qy = qs[i].qy;
		qs[i].q_idx = xy2d(n, qx, qy);
	}

	// sort the quadrants by their hilbert index
	std::sort(
		dpl::execution::par_unseq,
		qs.begin(),
		qs.end(),
		[](const Quad &a, const Quad &b) -> bool {
			return a.q_idx < b.q_idx;
		}
	);

}

/**
 * A microbenchmark measuring:
 * - runtime
 * - L1-3 cache misses
 * for the following workload:
 *
 * Given an input graph (represented using an edge list) with n vertices and m edges:
 * 1. permute the vertices of the graph using a random ordering
 * 	- this ensures that the edges will be distributed evenly across the adjacency matrix
 * 2. assign the edges of the graph to square quadrants.
 * 	- the following is a hparam that should be decided by the sizes of the L2, L3:
 * 		- side-length of a quadrant (0.25 * L2_SIZE) - must be a power of 2
 * 3. each quadrant contains the following data:
 * 	3.1 offset - the (x, y) coordinate to the top right corner of the quadrant
 * 		3.1.1 qx,
 * 		3.1.2 qy
 * 	3.2 nnz - the number of edges contained within the quadrant
 * 	3.3. data - a flattened array of the edges contained within the quadrant
 * 		e.g. [src_1, dest_1, src_2, dest_2, ..., src_nnz, dest_nnz]
 * 		The edges within a quadrant are ordered using the Hilbert Curve
 *
 * The benchmark iterates over the _quadrants_ in 2 ways:
 * I. Row major
 * II. Hilbert Order
 *
 * To ensure efficient sharing of L1, L2 caches among core threads, adjacent pairs of
 * quadrants are assigned to the pair of threads that belong to the same core.
 * e.g. given 6 cores with 2 threads per core and 16 quadrants:
-----------------------------------
HWThread        Thread        Core
0               0             0
1               0             1
2               0             2
3               0             3
4               0             4
5               0             5
6               1             0
7               1             1
8               1             2
9               1             3
10              1             4
11              1             5
 Threads will be assigned to quadrants as such:
 - Q1, Q2 -> HWThread 0, 6
 - Q3, Q4 -> HWThread 1, 7
 - ...
 - Q11, Q12 -> HWThread 5, 11
 *
 * This pairing of adjacent quadrants applies to both I. and II.
 *
 * Hypothesis:
 * Cache miss ratio, and as a result, execution time will be reduced for II. for graphs whose
 * vertex set (n) does not fit in memory
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
	opterr = 0;
	int opt;
	uint32_t n = 0; // n vertices
	uint64_t m = 0; // n edges
	std::string graph_name = "";
	std::string data_dir = "";
	bool hilbert = false;
	uint32_t q_side_len = 0;

	while ((opt = getopt(argc, argv, "hg:d:l:")) != -1) {
		switch (opt) {
			case 'h':
				hilbert = !hilbert;
				break;
			case 'g':
				graph_name = optarg;
				break;
			case 'd':
				data_dir = optarg;
				break;
			case 'l':
				q_side_len = atoi(optarg);
		}
	}
	topology_init();
	auto topo = get_cpuTopology();
	if (q_side_len > 0xFFFF || q_side_len == 0) {fmt::print("Invalid SideLength input: {}\n", q_side_len); return 0;}

	uint64_t L1_CACHE_SIZE = topo->cacheLevels[0].size;
	uint64_t L2_CACHE_SIZE = topo->cacheLevels[1].size;
	uint64_t L3_CACHE_SIZE = topo->cacheLevels[2].size;

//	fmt::print("L1_CACHE_SIZE: {}\n", L1_CACHE_SIZE);
//	fmt::print("L2_CACHE_SIZE: {}\n", L2_CACHE_SIZE);
//	fmt::print("L3_CACHE_SIZE: {}\n", L3_CACHE_SIZE);

	std::vector<uint32_t> iso_map;
	std::vector<std::pair<uint32_t, uint32_t>> edge_list;

	std::string graphs_dir = fmt::format("{}/{}", data_dir, "graphs");
	std::string graph_dir = fmt::format("{}/{}", graphs_dir, graph_name);
	std::string binary_order_path = fmt::format("{}/{}", graph_dir, "rnd.bin");
	std::string binary_edge_list_path = fmt::format("{}/{}", graph_dir, "comp.bin");
	read_binary_container<std::vector<uint32_t>>(binary_order_path, iso_map);
	read_binary_container<std::vector<std::pair<uint32_t, uint32_t>>>(binary_edge_list_path, edge_list);

	n = iso_map.size();
	m = edge_list.size();
//	fmt::print("n: {}\n", n);
//	fmt::print("m: {}\n", m);
	// remap the edges of the graph and sort by src, dest
#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < m; ++i) {
		uint32_t src = edge_list[i].first;
		uint32_t dest = edge_list[i].second;

		edge_list[i].first = iso_map[src];
		edge_list[i].second = iso_map[dest];
	}

	std::sort(dpl::execution::par_unseq, edge_list.begin(), edge_list.end());
	uint32_t n_quads_per_side = (n / q_side_len) + 1;
	uint64_t hyp_n = highest_power_of_2_greater_than(n);
	uint32_t hyp_n_quads = highest_power_of_2_greater_than(n_quads_per_side - 1);

	uint32_t n_vid_bits = log2(hyp_n);
	uint32_t n_q_bits = log2(hyp_n_quads);

//	fmt::print("q_side_len: {}\n", q_side_len);
//	fmt::print("n_quads_per_side: {}\n", n_quads_per_side);
//	fmt::print("n_vid_bits: {}\n", n_vid_bits);
//	fmt::print("n_q_bits: {}\n", n_q_bits);
	uint64_t total_n_quads = n_quads_per_side * n_quads_per_side;
	std::vector<Quad> qs(total_n_quads);
	std::vector<uint32_t> n_edges_seen_per_q(total_n_quads);

	// iterate over the edges of the graphs, computing the number of edges that will
	// be in each quadrant
	for (uint32_t i = 0; i < m; ++i) {
		uint32_t src = edge_list[i].first;
		uint32_t dest = edge_list[i].second;

		std::bitset<32> sbits = std::bitset<32>(src);
		std::bitset<32> dbits = std::bitset<32>(dest);

		// top n_q_bits of n_vid_bits of the vertex id gives the quadrant row and column given
		// source id and dest id, respectively
		unsigned q_mask, e_mask;
		short start_bit = n_vid_bits - n_q_bits;
		q_mask = ((1 << n_q_bits) - 1) << start_bit;
		e_mask = ((1 << start_bit) - 1);
		uint32_t qx = (src & q_mask) >> start_bit;
		uint32_t qy = (dest & q_mask) >> start_bit;
		uint32_t local_src = src & e_mask;
		uint32_t local_dest = dest & e_mask;

		uint32_t qidx = qx * n_quads_per_side + qy;
//		uint16_t t = 0;
		qs[qidx].qx = qx;
		qs[qidx].qy = qy;
		++qs[qidx].nnz;
	}
	uint64_t total_edges = 0;
#pragma omp parallel for schedule(static) reduction(+:total_edges)
	for (uint32_t i = 0; i < total_n_quads; ++i) {
		Quad &q = qs[i];
//		fmt::print("({} {}), {}\n", q.qx, q.qy, q.nnz);
		q.edges.resize(q.nnz * 2); // each non-zero consists of src, dest - so mult. by 2
		total_edges += q.edges.size();
	}

	assert(total_edges == m * 2);
	std::fill(n_edges_seen_per_q.begin(), n_edges_seen_per_q.end(), 0);

	// iterate over the edgelist (again) to populate the flat edge vector per quadrant
	for (uint32_t i = 0; i < m; ++i) {
		uint32_t src = edge_list[i].first;
		uint32_t dest = edge_list[i].second;

		unsigned q_mask, e_mask;
		short start_bit = n_vid_bits - n_q_bits;
		q_mask = ((1 << n_q_bits) - 1) << start_bit;
		e_mask = ((1 << start_bit) - 1);
		uint32_t qx = (src & q_mask) >> start_bit;
		uint32_t qy = (dest & q_mask) >> start_bit;
		uint32_t local_src = src & e_mask;
		uint32_t local_dest = dest & e_mask;
		uint32_t qidx = qx * n_quads_per_side + qy;
		uint32_t n_edges_seen = n_edges_seen_per_q[qidx];
		qs[qidx].edges[n_edges_seen] = local_src;
		qs[qidx].edges[n_edges_seen + 1] = local_dest;
		n_edges_seen_per_q[qidx] += 2;
	}
	std::string bin_file_name = hilbert ? "hqs" : "qs";
	std::string quad_path = fmt::format("{}/{}.bin", graph_dir, bin_file_name);

	if (hilbert) {
		// order the edge within a quadrant using the hilbert
		horder_edges_in_qs(qs, q_side_len);
		// order the quadrants using the hilbert order
		horder_qs(qs, n_quads_per_side);
	}

//	for (const auto &q: qs) {
//		fmt::print("q.qx, q.qy, q_idx, nnz: {} {} {} {}\n", q.qx, q.qy, q.q_idx, q.nnz);
//		fmt::print("q.edges: {}\n", q.edges);
//	}
	write_binary_qs(quad_path, qs, q_side_len, n_quads_per_side);

//#pragma omp parallel for schedule(static, 1)
//	for (uint32_t i = 0; i < 12; ++i) {
//		int cpu_num = sched_getcpu();
//		fmt::print("i: {} {} {}\n", i, omp_get_thread_num(), cpu_num);
//	}

	// prepare for pagerank
	// read the outdegree file
	std::string deg_path = fmt::format("{}/{}", graph_dir, "out_degs");
	std::vector<uint32_t> deg(n);
	std::vector<double> mapped_degs(n);
	read_text_degree_file(deg_path, deg);
#pragma omp parallel for
	for (uint32_t i = 0; i < n; ++i) {
		mapped_degs[iso_map[i]] = double(deg[i]);
	}
//	fmt::print("mapped_degs: {}\n", mapped_degs);

	double alpha = 0.85;
	const double init_score = 1.0f / n;
	const double base_score = (1.0f - alpha) / n;

	std::vector<double> outgoing_contrib(n);
	std::vector<double> scores(n);
	int num_iters = 20;

#pragma omp parallel for
	for (uint32_t i = 0; i < n; i++) {
		outgoing_contrib[i] = init_score / mapped_degs[i];
	}
	int n_threads = omp_get_max_threads();
	std::vector<pvector<double>> incoming_total(n_threads);  // private
	for (int t = 0; t < n_threads; ++t) {
		incoming_total[t].resize(n);
	}
//	fmt::print("n_threads: {}\n", n_threads);

	// init likwid markers
	LIKWID_MARKER_INIT;
	std::string likwid_marker_group_name = "Parallel-PageRank";
	const char * lgroup = likwid_marker_group_name.c_str();
#pragma omp parallel
	{
		LIKWID_MARKER_THREADINIT;
		LIKWID_MARKER_REGISTER(lgroup);
	}

//fmt::print("Iterating..\n");
	auto start_time = std::chrono::high_resolution_clock::now();

	for (int iter = 0; iter < num_iters; iter++) {
//		fmt::print("Iteration: {}\n", iter);
#pragma omp parallel
		{
			int tid = omp_get_thread_num();
			incoming_total[tid].fill(0);
//			fmt::print("incoming_total[tid]: {}\n", incoming_total[tid]);

LIKWID_MARKER_START(lgroup);
#pragma omp for schedule(static, 1)
			for (uint32_t quad_id = 0; quad_id < qs.size(); ++quad_id) {
				Quad &q = qs[quad_id];
				uint32_t qx = q.qx;
				uint32_t qy = q.qy;
				for (uint32_t i = 0; i < q.nnz; ++i) {
					uint32_t src = (qx * q_side_len) + q.edges[i * 2];
					uint32_t dest = (qy * q_side_len) + q.edges[i * 2 + 1];
					incoming_total[tid][dest] += outgoing_contrib[src];
				}
			}
LIKWID_MARKER_STOP(lgroup);

//#pragma omp barrier

#pragma omp for schedule(static)
			for (uint32_t i = 0; i < n; i++) {
				// read each threads incoming_total and sum
				double in_total = 0;
				for (int t = 0; t < n_threads; ++t) {
					in_total += incoming_total[t][i];
				}
				scores[i] = base_score + alpha * in_total;
				outgoing_contrib[i] = scores[i] / mapped_degs[i];
//				incoming_total[i] = 0;
			}
		}
	}
	LIKWID_MARKER_CLOSE;
	auto end_time = std::chrono::high_resolution_clock::now();
	uint64_t runtime = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

//	fmt::print("scores: {}\n", scores);
	std::string pr_path = fmt::format("{}/{}", graph_dir, "pr");
	std::vector<double> pr(n);
	std::vector<double> mapped_results(n);
#pragma omp parallel for
	for (uint32_t i = 0; i < n; ++i) {
		mapped_results[i] = scores[iso_map[i]];
	}
	read_binary_container<std::vector<double>>(pr_path, pr);
	bool valid_pr = std::equal(dpl::execution::par_unseq,
	                           mapped_results.begin(), mapped_results.end(),
	                           pr.begin(),
	                           [](double value1, double value2) {
		                           constexpr double epsilon = 1e-4;
		                           return std::fabs(value1 - value2) < epsilon;
	                           });
	fmt::print("valid_pr: {}\n", valid_pr);
//	for (uint32_t i = 0; i < 20; ++i) {
//		fmt::print("{} : {}\n", mapped_results[i], pr[i]);
//	}
	fmt::print("runtime: {}\n", runtime);
	return 0;
}