//
// Created by atrostan on 06/01/23.
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
#include <stdio.h>
#include "ips4o/ips4o.hpp"
#include "Quad.h"
#include <memory>

class PageRankGrid {
	uint8_t _n_threads;
	uint32_t _n_vertices;
//	std::unique_ptr<float[]> data;
	float *data = nullptr;

public:
	PageRankGrid(uint8_t n_threads, uint32_t n_vertices)
		: _n_threads{n_threads},
		  _n_vertices{n_vertices}
//		  data{std::make_unique<float[]>(n_threads * n_vertices)} {{
	{
		data = new float[n_threads * n_vertices]();
//		memset(data, 0, n_threads * n_vertices * sizeof(float));
	}

	uint8_t n_threads() const { return _n_threads; }

	uint32_t n_vertices() const { return _n_vertices; }

	float *operator[](uint8_t thread_id) { return data + thread_id * _n_vertices; }

	float &operator()(uint8_t thread_id, uint32_t vertex_id) {
		return data[thread_id * _n_vertices + vertex_id];
	}

	~PageRankGrid() {
		delete[]data;
	}
};

void merge_incoming_totals(PageRankGrid &grid, float *&incoming_totals, uint32_t start, uint32_t end,
                           uint8_t n_threads) {
#pragma omp parallel for
	for (uint32_t u = start; u < end; ++u)
		for (uint8_t t = 0; t < n_threads; ++t)
			incoming_totals[u] += grid(t, u);
}

void iterate_left_wing(Quad *qs, uint32_t n_qs, uint32_t wing_width, uint32_t q_side_len) {
#pragma omp parallel for schedule(static, 1)
	for (uint32_t i = 0; i < n_qs; ++i) {
		Quad &q = qs[i];
		for (uint32_t j = 0; j < q.nnz; ++j) {
			uint32_t src = q.edges[j * 2];
			uint32_t dest = q.edges[j * 2 + 1];
			uint32_t u = q.qx + src;
			uint32_t v = q.qy + dest;

		}
	}
}

/**
 * Given a graph name,
 * - Reads the parallel-slashburn reordered graph
 *   - stored in 3 binary files for the left wing, right wing, and tail
 *
 * - for each edge section (lw, rw, tail), iterates over the quadrants in parallel to compute the page rank
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
	opterr = 0;
	int opt;
	int num_expts = 0;
//	int n_threads = 0;
	std::string graph_name = "";
	std::string data_dir = "";
	uint32_t q_side_len = 0;
	uint32_t wing_width = 0;
	uint16_t num_iterations = 0;

	while ((opt = getopt(argc, argv, "g:d:l:e:w:i:")) != -1) {
		switch (opt) {
			case 'g':
				graph_name = optarg;
				break;
			case 'd':
				data_dir = optarg;
				break;
			case 'l':
				q_side_len = atoi(optarg);
			case 'w':
				wing_width = atoi(optarg);
			case 'i':
				num_iterations = atoi(optarg);
		}
	}
	assert(is_power_of_2(q_side_len));
	// requisite paths
	std::string graphs_dir = fmt::format("{}/{}", data_dir, "graphs");
	std::string graph_dir = fmt::format("{}/{}", graphs_dir, graph_name);
	std::string lw_path = fmt::format("{}/{}", graph_dir, "lw.bin");
	std::string rw_path = fmt::format("{}/{}", graph_dir, "rw.bin");
	std::string tail_path = fmt::format("{}/{}", graph_dir, "tail.bin");
	std::string out_deg_path = fmt::format("{}/{}", graph_dir, "out_degs");
	fmt::print("out_deg_path: {}\n", out_deg_path);

	// read the binary files containing the quadrants
	Quad *lw_qs = nullptr;
	Quad *rw_qs = nullptr;
	Quad *tail_qs = nullptr;
	std::vector<uint64_t> res;
	uint32_t n, n_lw_qs, n_rw_qs, n_tail_qs;
	uint64_t m;
	res = read_quad_array(lw_path, lw_qs, false);
	n_lw_qs = res[0];
	n = res[1];
	m = res[2];
	res = read_quad_array(rw_path, rw_qs, false);
	n_rw_qs = res[0];
	res = read_quad_array(tail_path, tail_qs, true);
	n_tail_qs = res[0];

	fmt::print("{} {} {} {} {}\n", n, m, n_lw_qs, n_rw_qs, n_tail_qs);

	int n_threads = omp_get_max_threads() / 2;
	fmt::print("n_threads: {}\n", n_threads);

	// read out degree file
	std::vector<uint32_t> out_degs(n);
	read_text_degree_file(out_deg_path, out_degs);

	std::vector<double> pr(n);
	std::vector<double> mapped_results(n);
	std::string pr_path = fmt::format("{}/{}", graph_dir, "pr");
	read_binary_container<std::vector<double>>(pr_path, pr);

	const float alpha = 0.85;
	const float init_score = 1.0f / n;
	const float base_score = (1.0f - alpha) / n;

	// thread private arrays
	PageRankGrid incoming_totals = PageRankGrid(n_threads, wing_width);

	// global arrays
	float *incoming_total = new float[n]();
	float *outgoing_contrib = new float[n]();
	float *scores = new float[n]();

	for (int iter = 0; iter < num_iterations; iter++) {

		// left wing
//		iterate_left_wing();

//		merge_incoming_totals()
		// right wing
//		iterate_right_wing();
		// tail
//		iterate_tail();

	}


	// clean up
	delete[]lw_qs;
	delete[]rw_qs;
	delete[]tail_qs;
}