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
#include "pvector.h"
#include "Quad.h"
#include "CSR.h"

struct vertex {
	ul id;
	ul degree;
};


template<class T, class I>
I lower_bound(T x, T *array, I l, I h) {
	// Invariant: This function will eventually return a value in      the range [begin, end]
	while (l < h) {
		int mid = l + (h - l) / 2;
		if (x <= array[mid]) {
			h = mid;
		} else {
			l = mid + 1;
		}
	}
	return l;
}

uint32_t get_min_id_in_range(uint32_t u, uint32_t v_min, uint32_t v_max, CSR &out_csr);

std::pair<uint32_t, uint32_t> edge_offset(uint32_t u, uint32_t v, uint32_t q_side_len);

uint32_t n_out_neighbours_between(uint32_t u, uint32_t v_min, uint32_t v_max, CSR &out_csr);

void horder_edges_in_q(Quad &q, uint32_t side_len);

uint32_t filter_empty_quads(Quad *&qs, uint32_t n_quads);

void normalize(std::vector<double> &orig, std::vector<double> &normalized);

bool sortByDescendingDegree(const vertex &lhs, const vertex &rhs);

void igraph_place_hubs(igraph_t &g, const int k, ul &hub_idx, std::vector<ul> &rank);

void create_isomorphism_map(std::vector<ul> &rank, std::map<ul, ul> &map);

void par_translate_edge_list(std::vector<Edge> &indexed_edges,
                             std::vector<Edge> &mapped_edges,
                             std::vector<ul> &iso_map, ull m);

bool is_power_of_2(uint32_t n);

void par_sort_edges(std::vector<Edge> &edges, Order ord, ul n);

uint32_t next_largest_multiple(uint32_t n, uint32_t critical_depth);

void print_seperator();

std::string dir_str(Direction d);

std::string corner_str(Corner c);

void print_quad(Quadrant &q);

std::pair<uint32_t, uint32_t> rotate_point(uint32_t cx, uint32_t cy, int angle, uint32_t x, uint32_t y);

uint32_t int_log(int base, uint32_t x);

template<typename T, typename U>
auto omp_accumulate(const pvector<T> &v, U init) {
	U sum = init;

#pragma omp parallel for reduction(+:sum)
	for (uint32_t i = 0; i < v.size(); i++)
		sum += v[i];

	return sum;
}

template<typename T>
void print_arr(T *arr, uint64_t size) {
	std::cout << "[";
	for (uint64_t i = 0; i < size - 1; ++i) { std::cout << arr[i] << ", "; }
	std::cout << arr[size - 1] << "]\n";
}

template<typename T>
void print_arr_slice(T *arr, uint64_t start, uint64_t end) {
	std::cout << "[";
	for (uint64_t i = start; i < end - 1; ++i) { std::cout << arr[i] << ", "; }
	std::cout << arr[end - 1] << "]\n";
}

#endif //GRAPH_PREPROCESS_UTIL_H
