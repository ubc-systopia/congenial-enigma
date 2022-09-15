//
// Created by atrostan on 30/08/22.
//

#ifndef GRAPH_PREPROCESS_UTIL_H
#define GRAPH_PREPROCESS_UTIL_H

#include "typedefs.h"
#include "igraph.h"
#include <map>
#include "Hilbert.h"

struct vertex {
	ul id;
	ul degree;
};

bool sortByDescendingDegree(const vertex &lhs, const vertex &rhs);

void igraph_place_hubs(igraph_t &g, const int k, ul &hub_idx, std::vector<ul> &rank);

void create_isomorphism_map(std::vector<ul> &rank, std::map<ul, ul> &map);

void par_translate_edge_list(std::vector<Edge> &indexed_edges,
                             std::vector<Edge> &mapped_edges,
                             std::vector<ul> &iso_map, ull m);

void par_sort_edges(std::vector<Edge> &edges, Order ord, ul n);

void print_quad(Quadrant &q);
#endif //GRAPH_PREPROCESS_UTIL_H
