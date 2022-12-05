//
// Created by atrostan on 29/11/22.
//

#include <getopt.h>
#include "stats.h"



/**
 * This executable will compute the algebraic stats of a directed graph using Eigen, Spectra
 * It expects an input graph directory to contain:
 *   - comp.bin - the edgelist of the graph (in binary format)
 *   - lcc.bin - the edgelist of the largest (weakly) connected component (LCC) of the graph
 *   - lscc.bin - the edgelist of the largest strongly connected component (LSCC) of the graph
 *
 * Given these edgelists, we compute:
 *
 * 1. Operator 2-Norm: the largest singular value of a directed graph
 * 2. Cyclic Eigenvalue: the largest absolute eigenvalue of the adjacency matrix of a directed graph
 * 3. Algebraic Connectivity: the second smallest eigenvalue of the Laplacian matrix of the matrix.
 *    Since the Laplacian can be defined for both the LCC and LSCC, we compute the Al. Conn. for both.
 *    For the directed LSCC, we use the Out-degree Laplacian.
 *    todo -> this is currently missing - can't verify correctness of smallest eval of lap mat
 * 4. Spectral Norm: the largest absolute eigenvalue of the graph's symmetric adjacency matrix.
 *    This requires symmetrizing the directed graph.
 * 5. Spectral Separation: the largest absolute eigenvalue of the adjacency matrix divided by the second largest
 *    absolute eigenvalue. This requires symmetrizing the directed graph.
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
	opterr = 0;
	int opt;
	uint32_t num_nodes = 0;
	std::string graph_dir = "";
	std::string sqlite3_db_path = "";

	while ((opt = getopt(argc, argv, "n:g:d:")) != -1) {
		switch (opt) {
			case 'n':
				num_nodes = atol(optarg);
				break;
			case 'g':
				graph_dir = optarg;
				break;
			case 'd':
				sqlite3_db_path = optarg;
				break;
		}
	}

	std::map<std::string, double> stats;

	const std::string edge_list_path = fmt::format("{}/comp.bin", graph_dir);
	const std::string lcc_edge_list_path = fmt::format("{}/lcc.bin", graph_dir);
	const std::string lscc_edge_list_path = fmt::format("{}/lscc.bin", graph_dir);

//	compute_eig_stats<double, std::vector>(false, true, edge_list_path, num_nodes, -1); // laplacian
//	compute_eig_stats<double, std::vector>(false, false, edge_list_path, num_nodes, 1);

	compute_eig_stats<double, std::vector>(false, false, edge_list_path, num_nodes, 1, sqlite3_db_path);

}