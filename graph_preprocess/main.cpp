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

	while ((opt = getopt(argc, argv, "g:b:m:o:")) != -1) {
		switch (opt) {
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
	std::vector<std::pair<ul, ul>> edges;
	std::vector<std::pair<ul, ul>> mapped_edges;
	edges.resize(num_edges);
	mapped_edges.resize(num_edges);

	boost::filesystem::path p(input_path);
	boost::filesystem::path dir = p.parent_path();
	std::string graph_name = dir.filename().string();

	std::pair<ul, ull> nm = read_edge_list(input_path, num_edges, edges, mapped_edges, graph_name, sqlite_db_path);
	ul n = nm.first;
	ull m = nm.second;
	fmt::print("n: {}\n", n);
	fmt::print("m: {}\n", m);
	std::ofstream outfile(output_path);

	// write the compressed, simplified edge list to file
	for (auto &kv: mapped_edges) {
		outfile << fmt::format("{} {}\n", kv.first, kv.second);
	}
	outfile.close();

	return 0;

}
