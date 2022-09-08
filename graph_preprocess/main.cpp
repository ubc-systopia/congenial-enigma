#include <iostream>
#include <getopt.h>
#include <string>
#include "io.h"
#include "typedefs.h"
#include <set>
#include <algorithm>    // std::random_shuffle
#include <cstdlib>      // std::rand, std::srand
#include <boost/filesystem.hpp>

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
				num_edges = atoi(optarg);
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
