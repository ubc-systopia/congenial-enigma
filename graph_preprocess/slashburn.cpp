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
#include <boost/unordered_set.hpp>
#include <unordered_set>
#include "order_slashburn.h"
#include "sql.h"
#include "command_line.h"
#include "benchmark.h"

int main(int argc, char *argv[]) {
	CLApp cli(argc, argv, "sequential slashburn");
	if (!cli.ParseArgs())
		return -1;
	Builder b(cli);

	Graph G = b.MakeGraph();
	float percent = cli.percent();

	ul num_vertices = G.num_nodes();
	ull num_edges = G.num_edges();
	bool directed = G.directed();
	std::string input_path = cli.filename();
	std::string output_path = cli.out_filename();
	std::string sqlite_db_path = cli.db_filename();

	/* turn on attribute handling */
	igraph_set_attribute_table(&igraph_cattribute_table);

	igraph_t g;
	igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
	igraph_add_vertices(&g, G.num_nodes(), NULL);
	igraph_vector_int_t flat_igraph_edges;

	uint64_t w = 0;
	boost::unordered_set<std::pair<uint32_t, uint32_t>> es;
	for (uint32_t u = 0; u < G.num_nodes(); ++u) {
		for (uint32_t v: G.out_neigh(u)) {
			uint32_t src, dest;
			if (u > v) { src = v; dest = u; }
			else { src = u; dest = v; }
			es.insert({src, dest});
		}
	}
	igraph_vector_int_init(&flat_igraph_edges, es.size() * 2);
	for (const auto &kv: es) {
			uint32_t u = kv.first;
			uint32_t v = kv.second;
			VECTOR(flat_igraph_edges)[w] = u;
			VECTOR(flat_igraph_edges)[w + 1] = v;
			w += 2;
	}
	fmt::print("\n");
	igraph_add_edges(&g, &flat_igraph_edges, NULL);
	igraph_vector_int_destroy(&flat_igraph_edges);
//	FILE *f;
//	f = fopen(input_path.c_str(), "r");

	// TODO slashburn impl ignores graph directedness
//	igraph_read_graph_edgelist(&g, f, num_vertices, 0);
//	std::vector<std::pair<ul, ul>> flat_edges(num_edges);
//	read_binary_edge_list_into_igraph(input_path, flat_edges, &g, num_vertices, num_edges, 0);
	ul n;
	ull i, m;

	n = igraph_vcount(&g);
	m = igraph_ecount(&g);
	igraph_vector_t v = IGRAPH_VECTOR_NULL;
	ul k = n * percent;
	if (k == 0) {
		k = 1;
	}
	fmt::print("g.directed: {}\n", g.directed);
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
		fmt::print("iter: {}, gcc.n: {}, hub_idx: {}, k: {}, spokes_end_idx: {}\n", iter, prev->n, hub_idx, k,
		           spokes_end_idx);
		res = order_igraph_slashburn(*prev, k, rank, hub_idx, spokes_end_idx);

		curr = get<0>(res);
		hub_idx = get<1>(res);
		spokes_end_idx = get<2>(res);

//		igraph_vit_t tvit;
//		igraph_vit_create(curr, igraph_vss_all(), &tvit);
//		while (!IGRAPH_VIT_END(tvit)) {
//			long int vid = (long int) IGRAPH_VIT_GET(tvit);
//			fmt::print("{} ", VAN(curr, "orig_id", vid));
//			IGRAPH_VIT_NEXT(tvit);
//		}
//		fmt::print("\n");


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
//		if (iter == 4) {
//			igraph_vit_t tvit;
//			igraph_vit_create(curr, igraph_vss_all(), &tvit);
//			while (!IGRAPH_VIT_END(tvit)) {
//				long int vid = (long int) IGRAPH_VIT_GET(tvit);
//				fmt::print("{}\n", VAN(curr, "orig_id", vid));
//				IGRAPH_VIT_NEXT(tvit);
//			}
//			fmt::print("\n");
//			break;
//		}

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

	single_val_set<int>(sqlite_db_path, "sb_k", "statistics", graph_name, k);
	single_val_set<int>(sqlite_db_path, "sb_n_iters", "statistics", graph_name, iter);
	single_val_set<int>(sqlite_db_path, "slashburn", "preproc", graph_name, int(slashburn_time.count()));

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