//
// Created by atrostan on 06/09/22.
//

#include <getopt.h>
#include "typedefs.h"
#include "io.h"

#include <chrono>
#include <fmt/core.h>

#include <fmt/ranges.h>

/* Only needed for the sake of this example. */
#include <iostream>
#include <thread>
#include <boost/config.hpp>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
#include <fmt/core.h>
#include <fstream>
#include <boost/filesystem.hpp>


int main(int argc, char *argv[]) {
	opterr = 0;
	int opt;
	ul num_nodes;
	ull num_edges;
	std::string input_path;
	std::string output_path;
	std::string sqlite_db_path;

	while ((opt = getopt(argc, argv, "n:m:p:g:b:o:")) != -1) {
		switch (opt) {
			case 'n':
				num_nodes = atoi(optarg);
				break;
			case 'm':
				num_edges = atoi(optarg);
				break;
			case 'g':
				input_path = optarg;
				break;
			case 'b':
				sqlite_db_path = optarg;
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

	using namespace boost;
	using namespace std;
	typedef adjacency_list<vecS, vecS, undirectedS,
			property<vertex_color_t, default_color_type,
					property<vertex_degree_t, int> > > Graph;
	typedef graph_traits<Graph>::vertex_descriptor Vertex;
	typedef graph_traits<Graph>::vertices_size_type size_type;

	boost::filesystem::path p(input_path);
	boost::filesystem::path dir = p.parent_path();
	std::string graph_name = dir.filename().string();


	std::string rev_cm_iso_path = fmt::format("{}/rev_cm", dir.string());
	std::string cm_iso_path = fmt::format("{}/cm", dir.string());


	std::vector<std::pair<ul, ul>> edges(num_edges);
	read_binary_edge_list(input_path, edges);

	Graph G(num_nodes);
	for (auto &kv: edges) {
		add_edge(kv.first, kv.second, G);
	}
	graph_traits<Graph>::vertex_iterator ui, ui_end;
	property_map<Graph, vertex_degree_t>::type deg = get(vertex_degree, G);
	for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
		deg[*ui] = degree(*ui, G);

	property_map<Graph, vertex_index_t>::type
			index_map = get(vertex_index, G);
//	std::cout << "original bandwidth: " << bandwidth(G) << std::endl;
	single_val_set_int(sqlite_db_path, "orig_bandwidth", "statistics", graph_name,  bandwidth(G));

	std::vector<Vertex> inv_perm(num_vertices(G));
	std::vector<size_type> perm(num_vertices(G));

	auto t1 = high_resolution_clock::now();
	cuthill_mckee_ordering(G, inv_perm.begin(), get(vertex_color, G),
	                       make_degree_map(G));
	auto t2 = high_resolution_clock::now();
	auto cm_runtime = duration_cast<time_unit>(t2 - t1);
	single_val_set_int(sqlite_db_path, "cuthill_mckee", "preproc", graph_name, int(cm_runtime.count()));
	ul cm_id = 0;
	ul rev_cm_id = 0;
	std::map<ul, ul> cm_map;
	std::map<ul, ul> rev_cm_map;

	// reverse cuthill mckee order
	for (std::vector<Vertex>::const_reverse_iterator i = inv_perm.rbegin();
	     i != inv_perm.rend(); ++i) {
		rev_cm_map[index_map[*i]] = rev_cm_id;
		rev_cm_id++;
	}

	// cuthill mckee order
	for (std::vector<Vertex>::const_iterator i = inv_perm.begin();
	     i != inv_perm.end(); ++i) {
		cm_map[index_map[*i]] = cm_id;
		cm_id++;
	}

	write_permutation(cm_iso_path, cm_map, num_nodes, num_edges);
	write_permutation(rev_cm_iso_path, rev_cm_map, num_nodes, num_edges);

	for (size_type c = 0; c != inv_perm.size(); ++c)
		perm[index_map[inv_perm[c]]] = c;
//	std::cout << "  bandwidth: "
//	          << bandwidth(G, make_iterator_property_map(&perm[0], index_map, perm[0]))
//	          << std::endl;

	single_val_set_int(sqlite_db_path, "cm_bandwidth", "statistics", graph_name,  bandwidth(G, make_iterator_property_map(&perm[0], index_map, perm[0])));

//	write_perm
}