#include <iostream>
#include <getopt.h>
#include <string>
#include "io.h"
#include "typedefs.h"
#include <set>
#include <algorithm>    // std::random_shuffle
#include <cstdlib>      // std::rand, std::srand
#include <boost/filesystem.hpp>
#include <eigen3/Eigen/Sparse>
#include "pvector.h"
#include <eigen3/Eigen/Core>

void convert_to_sparse_mat_and_save(bool laplacian, bool in_degree,
                                    std::vector<std::pair<uint32_t, uint32_t>> &edges,
                                    uint32_t n, std::string out_path) {
	typedef Eigen::Triplet<ul> T;
	pvector<uint32_t> deg(n, 0);
	uint64_t m = edges.size();
	uint64_t n_entries_in_mat = laplacian ? n + m : m;
	pvector<T> triples(n_entries_in_mat);
	uint32_t mat_val = laplacian ? -1 : 1;

#pragma omp parallel for schedule(static)
	for (uint64_t i = 0; i < m; ++i) {
		uint32_t src = edges[i].first;
		uint32_t dest = edges[i].second;
		triples[i] = T(src, dest, mat_val);
		if (in_degree) { deg[dest]++; }
		else { deg[src]++; }
	}
	if (laplacian) {
#pragma omp parallel for schedule(static)
		for (uint32_t i = 0; i < deg.size(); ++i) {
			triples[i + m] = T(i, i, deg[i]);
		}
	}

	Eigen::SparseMatrix<ul> M(n, n);
	M.setFromTriplets(triples.begin(), triples.end());
	write_binary_sparse(out_path, M);
}

/**
 * Read an undirected/directed graph from an edge-list file
 * Simplify the graph (Removes loop and/or multiple edges from the graph)
 * Compress the graph's ID space
 * @return
 */

int main(int argc, char *argv[]) {
	opterr = 0;
	int opt;
	ull num_edges = 0;
	std::string input_path;
	std::string sqlite_db_path;
	std::string output_path;
	std::vector<io_mode> io_modes;

	while ((opt = getopt(argc, argv, "tig:b:m:o:")) != -1) {
		switch (opt) {
			case 't':
				io_modes.push_back(text);
				break;
			case 'i':
				io_modes.push_back(binary);
				break;
			case 'g':
				input_path = optarg;
				break;
			case 'b':
				sqlite_db_path = optarg;
				break;
			case 'm':
				num_edges = atoll(optarg);
				break;
			case 'o':
				output_path = optarg;
				break;
			case '?':
				if (optopt == 'k')
					printf("Option -%c requires a long long.\n", optopt);
				else if (optopt == 'g' || optopt == 'o')
					printf("Option -%c requires a string.\n", optopt);
				else
					printf("Unknown option character '\\x%x'.\n", optopt);
				return 1;
			default:
				abort();
		}
	}

	// if no io mode supplied, default to text
	if (io_modes.empty()) {
		io_modes.push_back(text);
	}

	std::vector<std::pair<ul, ul>> mapped_edges;
	mapped_edges.resize(num_edges);

	boost::filesystem::path p(input_path);
	boost::filesystem::path dir = p.parent_path();
	std::string graph_name = dir.filename().string();

	std::vector<std::pair<ul, ul>> edges;
	edges.resize(num_edges);
	std::pair<ul, ull> nm = read_edge_list(input_path, num_edges, edges, mapped_edges, graph_name, sqlite_db_path);

//	std::pair<ul, ull> nm = par_read_edge_list(input_path, mapped_edges, graph_name, sqlite_db_path);
	ul n = nm.first;
	ull m = nm.second;
	fmt::print("n: {}\n", n);
	fmt::print("m: {}\n", m);

	for (auto &io_mode: io_modes) {
		write_edge_list(output_path, mapped_edges, io_mode);
	}
//
//	std::string sparse_mat_path = fmt::format("{}/{}", dir.string(), "mat.bin");
//	std::string laplacian_mat_path = fmt::format("{}/{}", dir.string(), "lap.bin");
//
//
//	// convert the compressed graph to a Sparse matrix and save
//	typedef Eigen::Triplet<ul> T;
//	pvector<T> triples(mapped_edges.size());
//	pvector<T> in_lap_triples(mapped_edges.size() + n);
//	pvector<T> out_lap_triples(mapped_edges.size() + n);
//	pvector<uint32_t> in_degs(n);
//	pvector<uint32_t> out_degs(n);
//
//#pragma omp parallel for schedule(static)
//	for (uint64_t i = 0; i < mapped_edges.size(); ++i) {
//		triples[i] = T(
//			mapped_edges[i].first, // src
//			mapped_edges[i].second, // dest
//			1 // weight
//		);
//	}
//	Eigen::SparseMatrix<ul> M(n, n);
//	M.setFromTriplets(triples.begin(), triples.end());

//	write_binary_sparse()
//	Eigen::SparseMatrix::
//

	// DEBUG
//	std::vector<std::pair<ul, ul>> tmp_edges(m);
//	read_binary_edge_list(
//			fmt::format("{}.bin", output_path),
//			tmp_edges
//	);
//	for (auto &e: tmp_edges) {
//		fmt::print("e: {}\n", e);
//	}

	return 0;

}
