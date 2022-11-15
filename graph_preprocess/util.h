//
// Created by atrostan on 30/08/22.
//

#ifndef GRAPH_PREPROCESS_UTIL_H
#define GRAPH_PREPROCESS_UTIL_H

#include "typedefs.h"
#include "igraph.h"
#include <map>
//#include "furhilbert.h"
#include "Hilbert.h"

struct vertex {
	ul id;
	ul degree;
};

void normalize(std::vector<double> &orig, std::vector<double> &normalized);

bool sortByDescendingDegree(const vertex &lhs, const vertex &rhs);

void igraph_place_hubs(igraph_t &g, const int k, ul &hub_idx, std::vector<ul> &rank);

void create_isomorphism_map(std::vector<ul> &rank, std::map<ul, ul> &map);

void par_translate_edge_list(std::vector<Edge> &indexed_edges,
                             std::vector<Edge> &mapped_edges,
                             std::vector<ul> &iso_map, ull m);

void par_sort_edges(std::vector<Edge> &edges, Order ord, ul n);

uint32_t next_largest_multiple(uint32_t n, uint32_t critical_depth);

void print_seperator();

std::string dir_str(Direction d);

std::string corner_str(Corner c);

void print_quad(Quadrant &q);

std::pair<uint32_t, uint32_t> rotate_point(uint32_t cx, uint32_t cy, int angle, uint32_t x, uint32_t y);

uint32_t int_log(int base, uint32_t x);

#endif //GRAPH_PREPROCESS_UTIL_H
