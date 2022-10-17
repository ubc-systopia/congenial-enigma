//
// Created by atrostan on 30/08/22.
//

#include <tuple>
#include "order_slashburn.h"

#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/numeric>
#include <iostream>

bool descSubgraphSize(const igraph_t *subgraph1, const igraph_t *subgraph2) {
	return subgraph1->n > subgraph2->n;
}

std::tuple<igraph_t *, ul, ul>
order_igraph_slashburn(igraph_t &g, const int k, std::vector<ul> &rank, ul &hub_idx, ul &spokes_end_idx) {
	// remove k hubs
	igraph_place_hubs(g, k, hub_idx, rank);

	ul num_components;

	igraph_graph_list_t complist;

	igraph_graph_list_init(&complist, 0);

	igraph_decompose(&g, &complist, IGRAPH_WEAK, -1, 0);

	num_components = igraph_graph_list_size(&complist);
	if (num_components == 0) {
		return std::tuple(&g, -1, -1);
	}

	// sort the component list by descending order of size to identify
	// the gcc and spokes
	std::vector<igraph_t *> subgraphs(num_components);
	ul i;

	for (i = 0; i < igraph_graph_list_size(&complist); ++i) {
		subgraphs[i] = igraph_graph_list_get_ptr(&complist, i);
		igraph_vit_t viter;
		igraph_vit_create(subgraphs[i], igraph_vss_all(), &viter);
	}

	std::sort(dpl::execution::par_unseq, subgraphs.begin(), subgraphs.end(), descSubgraphSize);

	igraph_t *gcc = subgraphs[0];

	std::vector<igraph_t *> spokes = std::vector<igraph_t *>(
			subgraphs.begin() + 1, subgraphs.end()
	);

	ul spokes_len = std::transform_reduce(
			dpl::execution::par_unseq,
			spokes.begin(), spokes.end(), 0,
			[](ul l, ul r) { return l + r; },
			[](igraph_t *c) { return c->n; }
	);

	ul spoke_idx = spokes_end_idx - spokes_len;
	for (auto spoke: spokes) {
		// iterate through the vertices in each spoke and place them in the end of the permutation array
		igraph_vit_t vit;
		igraph_vit_create(spoke, igraph_vss_all(), &vit);
//		std::cout<< "break\n";
		while (!IGRAPH_VIT_END(vit)) {

			long int vid = (long int) IGRAPH_VIT_GET(vit);
//			std::cout << "vid: " << VAN(spoke, "orig_id", vid) << "\n";
			rank[spoke_idx] = VAN(spoke, "orig_id", vid);
			IGRAPH_VIT_NEXT(vit);
			++spoke_idx;
		}

	}

	// all spoke vertices have been written to the permutation array, so destroy the ptr vectors (still need the gcc ptr, though)
	for (i = 1; i < igraph_graph_list_size(&complist); ++i) {
		auto s1 = subgraphs[i];
		igraph_destroy(s1);
//		igraph_free(s1); //todo apparently no longer needed..
	}

//	igraph_graph_list_destroy(&complist);

	auto res = std::make_tuple(gcc, hub_idx, spokes_end_idx - spokes_len);
	return res;
}