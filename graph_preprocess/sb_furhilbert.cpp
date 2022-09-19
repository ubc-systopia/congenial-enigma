//
// Created by atrostan on 14/09/22.
//


/**
* given a graph that has been reordered using slashburn, compute
 * the fur-hilbert order of the edges contained within it's k*i core
 * (upper, left square in the adjacency matrix)
*/

#include <random>
#include "sb_furhilbert.h"

int main(int argc, char *argv[]) {
	opterr = 0;
	int opt;
	uint32_t n = 0;
	uint32_t m = 0;
	while ((opt = getopt(argc, argv, "n:")) != -1) {
		switch (opt) {
			case 'n':
				n = atoi(optarg);
				break;
//			case 'm':
//				m = atoi(optarg);
//				break;
		}
	}
	uint32_t imin;
	uint32_t imax;
	uint32_t jmin;
	uint32_t jmax;
	uint32_t i, j;
	uint64_t hilbert_counter = 0;
	i = j = 0;

	/**
	 * given the number of available threads, compute the next largest value of n
	 * that can be evenly split among the threads
   * e.g. given 20 threads, the adjacency matrix will be split into 32 quadrants,
   * (Critical depth = 5), to ensure that all threads are busy computing the furhilbert
   * curve.
   * If n = 215 and critical_depth = 5, the size of the square we should recusively fill
   * with the furhilbert curve is next_largest_multiple(215, 5) = 224
   * This rounding produces a constant amount of additional work but simplifies
   * the furhilbert parallel traversal logic
	 */

	int n_threads = omp_get_max_threads();
//	int n_threads = m;
//	int n_threads = 5;
	int n_quads = 1;
	int critical_depth = 0;
	while (n_quads < n_threads) {
		critical_depth += 1;
		n_quads = n_quads * 4;
	}
	fmt::print("critical_depth: {}\n", critical_depth);
	n = next_largest_multiple(n, critical_depth + 2);
	fmt::print("n: {}\n", n);
	fmt::print("n_threads: {}\n", n_threads);
	fmt::print("n_quads: {}\n", n_quads);

	std::vector<Quadrant> qs(n_quads);

	QuadrantCalculator qc = QuadrantCalculator(critical_depth, n, qs);

	std::string debug_path = fmt::format(
			"/home/atrostan/Workspace/repos/congenial-enigma/graph_preprocess/cmake-build-debug/debug/hilbert"
	);

	std::ofstream outfile(debug_path);
//	qs[0].rot = 90;
	for (auto &q: qs) {
		int i = q.start_x;
		int j = q.start_y;
		imin = q.start_x;
		jmin = q.start_y;
		imax = q.end_x;
		jmax = q.end_y;
		uint32_t mid_i = ((q.start_x + q.end_x) / 2) ;
		uint32_t mid_j = ((q.start_y + q.end_y) / 2) ;
		print_quad(q);
//		fmt::print("mid_i: {}, mid_j: {}\n", mid_i, mid_j);
		FUR_HILBERT_FOR(i, j, imin, imax, jmin, jmax)
							{
								std::pair<uint32_t, uint32_t> p = rotate_point(mid_j, mid_i, q.rot, j, i);
//								std::pair<uint32_t, uint32_t> p = std::make_pair(i, j);
								outfile << p.second << " " << p.first << "\n";
							}
		FUR_HILBERT_END(i, j);

	}
	outfile.close();

//
//	// test
//
//
//	std::random_device rd; // obtain a random number from hardware
//	std::mt19937 n_threads_gen(rd()); // seed the generator
//	std::mt19937 n_gen(rd()); // seed the generator
//	std::uniform_int_distribution<> n_threads_distr(4, 129); // define the range
//	std::uniform_int_distribution<> n_distr(50, 100'000); // define the range
//
//

//	for (int m = 0; m < 40; ++m) {
//		print_seperator();
//		int n_threads = n_threads_distr(n_threads_gen);
//		n = n_distr(n_gen);
//		fmt::print("n_threads: {}\n", n_threads);
//		fmt::print("n: {}\n", n);
//		int n_quads = 1;
//		int critical_depth = 0;
//		while (n_quads < n_threads) {
//			critical_depth += 1;
//			n_quads = n_quads * 4;
//		}
//		n = next_largest_multiple(n, critical_depth + 2);
//		std::vector<Quadrant> qs(n_quads);
//
//		QuadrantCalculator qc = QuadrantCalculator(critical_depth, n, qs);
//		for (auto &q: qs) {
//			print_quad(q);
//		}
//		print_seperator();
//	}
}


