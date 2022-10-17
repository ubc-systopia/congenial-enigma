//
// Created by atrostan on 30/08/22.
//
#include <getopt.h>
#include "typedefs.h"
#include "io.h"
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>    // std::random_shuffle
#include <cstdlib>      // std::rand, std::srand
#include <cassert>
#include <igraph.h>
#include <boost/filesystem/path.hpp>
#include "order_slashburn.h"

int main(int argc, char *argv[]) {
	opterr = 0;
	int opt;
	ul num_vertices;
	ull num_edges;
	float percent;
	bool directed = false;
	std::string input_path;
	std::string output_path;
	std::string sqlite_db_path;

	while ((opt = getopt(argc, argv, "dn:m:p:g:b:o:")) != -1) {
		switch (opt) {
			case 'd':
				directed = !directed;
				break;
			case 'n':
				num_vertices = atol(optarg);
				break;
			case 'm':
				num_edges = atoll(optarg);
				break;
			case 'p':
				percent = atof(optarg);
				break;
			case 'g':
				input_path = optarg;
				break;
			case 'b':
				sqlite_db_path = optarg;
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
	/* turn on attribute handling */
	igraph_set_attribute_table(&igraph_cattribute_table);


	igraph_t g;
//	FILE *f;
//	f = fopen(input_path.c_str(), "r");

	// TODO slashburn impl ignores graph directedness
//	igraph_read_graph_edgelist(&g, f, num_vertices, 0);
	std::vector<std::pair<ul, ul>> flat_edges(num_edges);
	read_binary_edge_list_into_igraph(input_path, flat_edges, &g, num_vertices, num_edges, 0);
	ul n;
	ull i, m;

	n = igraph_vcount(&g);
	m = igraph_ecount(&g);
	igraph_vector_t v = IGRAPH_VECTOR_NULL;
	ul k = n * percent;
	k = 1; // tmp
	if (k == 0) {
		k = 1;
	}
	fmt::print("k: {}\n", k);
	fmt::print("n, m: {} {} \n", n, m);

	// set a boolean and numeric attribute for each vertex
	// the original id - this needs to be maintained as we modify the graph
	// seen - whether we've already assigned this vertex a final rank id
	for (i = 0; i < n; i++) {
		igraph_vs_t vs;
		igraph_vit_t vit;
		igraph_vs_adj(&vs, i, IGRAPH_ALL);
		igraph_vit_create(&g, vs, &vit);
		SETVAN(&g, "orig_id", i, i);
		SETVAB(&g, "seen", i, false);
	}
	vector<ul> rank;
	rank.resize(n, -1);
	ul hub_idx = 0;
	ul spokes_end_idx = n;
	// iterative result of slashburn; a triple of  (n_nodes in gcc, hub idx, spokes end idx)
	std::tuple<igraph_t *, ul, ul> res;
	igraph_t *prev;
	igraph_t *curr;

	prev = &g;
	ul iter = 0;
	auto start = std::chrono::high_resolution_clock::now();
	while (prev->n >= k) {
		fmt::print("gcc.n: {}, hub_idx: {}, k: {}, spokes_end_idx: {}\n", prev->n, hub_idx, k, spokes_end_idx);
		res = order_igraph_slashburn(*prev, k, rank, hub_idx, spokes_end_idx);

		curr = get<0>(res);
		hub_idx = get<1>(res);
		spokes_end_idx = get<2>(res);

		// if the greatest connected component's size is lesser than k, place the vertices in the gcc
		// after the latest placed hubs
		if (curr->n < k) {
			igraph_vit_t vit;
			igraph_vit_create(curr, igraph_vss_all(), &vit);
			while (!IGRAPH_VIT_END(vit)) {
				long int vid = (long int) IGRAPH_VIT_GET(vit);
				rank[hub_idx] = VAN(curr, "orig_id", vid);
				IGRAPH_VIT_NEXT(vit);
				++hub_idx;
			}
		}

		// destroy the previous gcc and reassign the pointers
		igraph_destroy(prev);
		prev = curr;
		iter += 1;
	}

	double wing_width_ratio = double(k * iter) / n;
	boost::filesystem::path p(input_path);
	boost::filesystem::path dir = p.parent_path();
	std::string graph_name = dir.filename().string();

	auto end = std::chrono::high_resolution_clock::now();
	auto slashburn_time = duration_cast<time_unit>(end - start);

	single_val_set_int(sqlite_db_path, "sb_k", "statistics", graph_name, k);
	single_val_set_int(sqlite_db_path, "sb_n_iters", "statistics", graph_name, iter);
	single_val_set_int(sqlite_db_path, "slashburn", "preproc", graph_name, int(slashburn_time.count()));

	fmt::print("rank: {}\n", rank);
	// all vertices have been assigned an index in the slashburn ordering
	assert(std::none_of(rank.begin(), rank.end(), [](ul v) {
		return v == std::numeric_limits<ul>::max();
	}));

	std::map<ul, ul> iso_map;
	create_isomorphism_map(rank, iso_map);

	write_permutation(output_path, iso_map, n, m);

	// DESTROYS!!
	igraph_destroy(&g);
	return 0;

}