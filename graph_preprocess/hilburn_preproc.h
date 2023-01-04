//
// Created by atrostan on 28/12/22.
//
#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/numeric>

#include <getopt.h>
#include <cstdint>
#include <vector>
#include <string>
#include <omp.h>
#include <likwid-marker.h>
#include <likwid.h>
#include "fmt/core.h"
#include "fmt/ranges.h"
#include "io.h"
#include "omp.h"
#include "util.h"
#include <bitset>
#include <limits.h>
#include "ips4o/ips4o.hpp"

#ifndef GRAPH_PREPROCESS_HILBURN_PREPROC_H
#define GRAPH_PREPROCESS_HILBURN_PREPROC_H

typedef uint32_t VertexID;
typedef uint64_t EdgeID;
typedef uint64_t QuadID;
typedef float ScoreT;

class Quad {
public:
	uint32_t qx;
	uint32_t qy;
	uint32_t q_idx;
	uint32_t nnz;
	uint32_t *edges = nullptr;

	Quad() {
		nnz = 0;
	}

	Quad(Quad &&other) : qx(other.qx), qy(other.qy),
	                     q_idx(other.q_idx), nnz(other.nnz), edges(other.edges) {
		other.qx = -1;
		other.qy = -1;
		other.q_idx = -1;
		other.nnz = -1;
		other.edges = nullptr;
	}

	Quad &operator=(Quad &&other) {
		if (this != &other) {
			ReleaseResources();
			qx = other.qx;
			qy = other.qy;
			q_idx = other.q_idx;
			nnz = other.nnz;
			edges = other.edges;
			other.qx = -1;
			other.qy = -1;
			other.q_idx = -1;
			other.nnz = -1;
			other.edges = nullptr;
		}
		return *this;
	}

	void ReleaseResources() {
		if (edges != nullptr) {
			delete[] edges;
		}
	}

	~Quad() {
		ReleaseResources();
	}
};

class CSR {
public:
	uint64_t *index = nullptr;
	uint32_t *neighbours = nullptr;
	uint32_t num_nodes;
	uint64_t num_edges;
	bool in;

	CSR(uint32_t n, uint64_t m, bool i) {
		in = i;
		num_nodes = n;
		num_edges = m;
		index = new uint64_t[num_nodes + 1]();
		neighbours = new uint32_t[num_edges]();
	}

	~CSR() {
		if (index != nullptr) { delete[] index; }
		if (neighbours != nullptr) { delete[] neighbours; }
	}

	void print_contents() {
		fmt::print("Index: [");
		for (uint32_t i = 0; i < num_nodes; ++i) { fmt::print("{}, ", index[i]); }
		fmt::print("]\n");
		fmt::print("Neighbours: [");
		for (uint32_t i = 0; i < num_edges; ++i) { fmt::print("{}, ", neighbours[i]); }
		fmt::print("]\n");
	}

	/**
	 * Input edge list assumed to be sorted by either:
	 * (src, dest) -> for out-csr
	 * (dest, src) -> for in-csr
	 * Correspondingly, input degrees array assumed to be either: out/in degrees.
	 * @param degrees
	 * @param edge_list
	 */
	void par_populate(std::vector<uint32_t> &degrees, std::vector<std::pair<uint32_t, uint32_t>> &edge_list) {
		uint64_t offset = 0;

		std::exclusive_scan(
			dpl::execution::par_unseq,
			degrees.begin(),
			degrees.end(),
			index,
			0,
			[](uint32_t l, uint32_t r) -> uint32_t { return l + r; }
		);

#pragma omp parallel for schedule(static)
		for (uint32_t u = 0; u < num_nodes; ++u) {
			uint32_t degree = degrees[u];
			uint32_t neigh_count = 0;
			for (uint32_t i = 0; i < degree; ++i) {
				uint32_t v = -1;
				if (in) { v = edge_list[index[u] + i].first; }
				else { v = edge_list[index[u] + i].second; }
				neighbours[index[u] + neigh_count] = v;
				++neigh_count;
			}
		}
		index[num_nodes] = num_edges;
	}

};


void read_degs(std::string out_path, std::string in_path, std::vector<uint32_t> &out_degs,
               std::vector<uint32_t> &in_degs, uint32_t n, std::vector<uint32_t> &iso_map);

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

std::pair<uint64_t, uint64_t> binary_search(uint32_t *arr, uint64_t x, uint64_t l, uint64_t h) {
	uint64_t mid;
	while (h - l > 1) {
		mid = (h + l) / 2;
		if (arr[mid] < x) l = mid + 1;
		else h = mid;
	}
	return {l, h};
//	if (arr[l] == x) return l;
//	else if (arr[h] == x) return h;
//	else return -1;
}


#endif //GRAPH_PREPROCESS_HILBURN_PREPROC_H
