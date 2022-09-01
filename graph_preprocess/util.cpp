//
// Created by atrostan on 30/08/22.
//

#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include <vector>
#include <map>
#include <fstream>
#include "util.h"


bool sortByDescendingDegree(const vertex &lhs, const vertex &rhs) {
	if (lhs.degree != rhs.degree)
		return lhs.degree > rhs.degree;
	return lhs.id < rhs.id;
}

void igraph_place_hubs(igraph_t &g, const int k, ul &hub_idx, std::vector<ul> &rank) {
	igraph_vector_t igraph_deg;
	igraph_vector_init(&igraph_deg, 0);
	ul i;

	// select k highest degree vertices
	igraph_degree(&g, &igraph_deg, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
	// use basic hub ordering for now TODO
	std::vector<vertex> vertices(g.n);
	for (ul j = 0; j < g.n; ++j) {
		vertex v;
		v.id = j;
		v.degree = VECTOR(igraph_deg)[j];
		vertices[j] = v;
	}

	std::sort(dpl::execution::par_unseq, vertices.begin(), vertices.end(), sortByDescendingDegree);

	igraph_vector_t hubs_to_remove;
	igraph_vector_init(&hubs_to_remove, k);
	ul j;
	for (j = 0; j < k; ++j) {
		vertex v = vertices[j];
		VECTOR(hubs_to_remove)[j] = v.id;
		rank[hub_idx + j] = VAN(&g, "orig_id", v.id);
	}
	hub_idx += k;

	igraph_delete_vertices(&g, igraph_vss_vector(&hubs_to_remove));

	igraph_vector_destroy(&igraph_deg);
	igraph_vector_destroy(&hubs_to_remove);
}



void create_isomorphism_map(std::vector<ul> &rank, std::map<ul, ul> &map) {
	for (ul i = 0; i < rank.size(); ++i) {
		map[rank[i]] = i;
	}
}

