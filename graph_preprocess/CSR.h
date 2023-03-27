//
// Created by atrostan on 11/01/23.
//

#include "util.h"
#include "fmt/core.h"
#include "fmt/ranges.h"
#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/numeric>
#ifndef GRAPH_PREPROCESS_CSR_H
#define GRAPH_PREPROCESS_CSR_H



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


#endif //GRAPH_PREPROCESS_CSR_H
