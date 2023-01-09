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
	uint16_t _n_threads;
	uint32_t _n_vertices;
//	std::unique_ptr<float[]> data;
	float *data = nullptr;

public:
	PageRankGrid(uint16_t n_threads, uint32_t n_vertices)
		: _n_threads{n_threads},
		  _n_vertices{n_vertices}
//		  data{std::make_unique<float[]>(n_threads * n_vertices)} {{
	{
		// todo possibly look into omp_alloc to ensure exclusive access to regions
		// by cores
		data = new float[n_threads * n_vertices]();
//		memset(data, 0, n_threads * n_vertices * sizeof(float));
	}

	uint16_t n_threads() const { return _n_threads; }

	uint32_t n_vertices() const { return _n_vertices; }

	float *operator[](uint16_t thread_id) { return data + thread_id * _n_vertices; }

	float &operator()(uint16_t thread_id, uint32_t vertex_id) {
//		return data[thread_id * _n_vertices + vertex_id];
		return data[vertex_id * _n_threads + thread_id];
	}

	~PageRankGrid() {
		if (data != nullptr) {
			delete[]data;
		}

	}
};

void merge_incoming_totals(PageRankGrid &grid, float *incoming_total,
                           uint32_t offset, uint32_t end,
                           uint16_t n_threads) {
// todo simd
//#pragma omp parallel for simd schedule(static)
#pragma omp parallel for schedule(static)
	for (uint32_t u = 0; u < end; ++u)
		for (uint16_t t = 0; t < n_threads; ++t) {
//			fmt::print("offset + u: {}\n", offset + u);
//			fmt::print("grid(t, u): {}\n", grid(t, u));
			incoming_total[offset + u] += grid(t, u);
			grid(t, u) = 0;
		}


}

void iterate_left_wing(Quad *qs, uint32_t n_qs,
                       uint32_t q_side_len,
                       PageRankGrid &grid, float *outgoing_contrib,
                       bool debug) {
	std::vector<uint16_t> thread_assignments;
	// if debugging, keep track of which thread computed which quadrant
	if (debug) {
		thread_assignments.resize(n_qs);
	}
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
//#pragma omp for schedule(runtime)
#pragma omp for schedule(static, 1)
//#pragma omp for schedule(dynamic, 1)
		for (uint32_t i = 0; i < n_qs; ++i) {
			Quad &q = qs[i];
			if (debug) {
				thread_assignments[i] = tid;
			}
			//todo would the following benefit from vectorization? (test)
//#pragma omp simd
			for (uint32_t j = 0; j < q.nnz; ++j) {
				uint32_t src = q.edges[j * 2];
				uint32_t dest = q.edges[j * 2 + 1];
				uint32_t u = q_side_len * q.qx + src;
				uint32_t v = q_side_len * q.qy + dest;
//				grid[tid][v] += outgoing_contrib[u];
				grid(tid, v) += outgoing_contrib[u];

			}
		}
	}
	if (debug) {
		fmt::print("thread_assignments: {}\n", thread_assignments);
		omp_sched_t sched;
		int chunk;
		omp_get_schedule(&sched, &chunk);
		fmt::print("sched, chunk: {} {}\n", sched, chunk);
	}


}

void iterate_tail(Quad *qs, uint32_t n_qs,
                  uint32_t q_side_len,
                  float *incoming_total, float *outgoing_contrib,
                  uint32_t n) {
	// each rectangle in the tail is guaranteed to write to an exclusive
	// range of destination ids
	// rectangles are sorted by qy
#pragma omp parallel for schedule(static, 1)
	for (uint32_t i = 0; i < n_qs; ++i) {
		Quad &q = qs[i];
		uint32_t v_end = q.qy + q_side_len > n ? n : q.qy + q_side_len;
		for (uint32_t j = 0; j < q.nnz; ++j) {
			uint32_t src = q.edges[j * 2];
			uint32_t dest = q.edges[j * 2 + 1];

//					omp_all
			uint32_t u = q.qx + src;
			// no need to offset destination ids - offset included in q.qy
			uint32_t v = q.qy + dest;

			// safe to directly write to global array
			incoming_total[v] += outgoing_contrib[u];
		}
	}

	return;
}

void iterate_right_wing(Quad *qs, uint32_t n_qs,
                        uint32_t wing_width, uint32_t q_side_len,
                        PageRankGrid &grid, float *outgoing_contrib,
                        uint32_t *cumulative_n_qs_per_right_wing_stripe,
                        uint32_t n_stripes_in_right_wing,
                        uint32_t n_qs_per_stripe,
                        float *incoming_total,
                        uint32_t n, uint16_t n_threads) {
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		for (uint32_t i = 0; i < n_stripes_in_right_wing; ++i) {
			uint32_t stripe_start = cumulative_n_qs_per_right_wing_stripe[i];
			uint32_t stripe_end = cumulative_n_qs_per_right_wing_stripe[i + 1];
#pragma omp for schedule(static, 1)
			for (uint32_t j = stripe_start; j < stripe_end; ++j) {
				Quad &q = qs[j];
//				fmt::print("q.nnz: {}\n", q.nnz);
				// the quad column should be zero-based to correctly index into thread
				// local array
				// todo maybe this should be done during preprocessing
				uint32_t qy = q.qy % n_qs_per_stripe;
				for (uint32_t k = 0; k < q.nnz; ++k) {
					uint32_t src = q.edges[k * 2];
					uint32_t dest = q.edges[k * 2 + 1];

//					omp_all
					uint32_t u = q_side_len * q.qx + src;
					uint32_t v = q_side_len * qy + dest;
					grid(tid, v) += outgoing_contrib[u];
//					grid[tid][v] += outgoing_contrib[u];
				}
			}

			// once we've completed this stripe, merge the results into
			// the global incoming total array
			uint32_t v_offset = wing_width + (n_qs_per_stripe * q_side_len * i);
			uint32_t v_end = wing_width + (n_qs_per_stripe * q_side_len * (i + 1));
			v_end = (v_end > n) ? n : v_end;
			v_end -= v_offset;
			merge_incoming_totals(grid, incoming_total, v_offset, v_end, n_threads);

//			fmt::print("{} {} {} {}\n", stripe_start, stripe_end, v_offset, v_end);
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
	bool debug = false;
	std::string schedule = ""; // 1 of {"static", "dynamic"}

	while ((opt = getopt(argc, argv, "bg:d:l:e:w:i:s:")) != -1) {
		switch (opt) {
			case 'b':
				debug = !debug;
				break;
			case 'g':
				graph_name = optarg;
				break;
			case 'd':
				data_dir = optarg;
				break;
			case 'l':
				q_side_len = atoi(optarg);
				break;
			case 'w':
				wing_width = atoi(optarg);
				break;
			case 'i':
				num_iterations = atoi(optarg);
				break;
			case 's':
				schedule = optarg;
				break;
		}
	}
	assert(is_power_of_2(q_side_len));

//	if (schedule == "") {
//		std::cout << "Please supply a valid OMP Schedule string (1 of: {static, dynamic}\n";
//		return 1;
//	}
//	if (schedule == "dynamic") { omp_set_schedule(omp_sched_dynamic, 1); }
//	else if (schedule == "static") { omp_set_schedule(omp_sched_static, 1); }
//	else { omp_set_schedule(omp_sched_static, 1); } //default

	// requisite input paths
	std::string graphs_dir = fmt::format("{}/{}", data_dir, "graphs");
	std::string graph_dir = fmt::format("{}/{}", graphs_dir, graph_name);
	std::string lw_path = fmt::format("{}/{}", graph_dir, "lw.bin");
	std::string rw_path = fmt::format("{}/{}", graph_dir, "rw.bin");
	std::string tail_path = fmt::format("{}/{}", graph_dir, "tail.bin");
	std::string out_deg_path = fmt::format("{}/{}", graph_dir, "out_degs");
	std::string binary_order_path = fmt::format("{}/{}", graph_dir, "parsb.bin");

	std::vector<uint32_t> iso_map;
	read_binary_container<std::vector<uint32_t>>(binary_order_path, iso_map);

	fmt::print("out_deg_path: {}\n", out_deg_path);

	// read the binary files containing the quadrants
	Quad *lw_qs = nullptr;
	Quad *rw_qs = nullptr;
	Quad *tail_qs = nullptr;
	uint32_t *cumulative_n_qs_per_right_wing_stripe;
	uint32_t *empty;
	std::vector<uint64_t> res;
	uint32_t n, n_lw_qs, n_rw_qs, n_tail_qs;
	uint64_t m;
	res = read_quad_array(lw_path, lw_qs, false, false, empty);
	n_lw_qs = res[0];
	n = res[1];
	m = res[2];
	res = read_quad_array(rw_path, rw_qs, false, true, cumulative_n_qs_per_right_wing_stripe);
	n_rw_qs = res[0];
	res = read_quad_array(tail_path, tail_qs, true, false, empty);
	n_tail_qs = res[0];

	if (debug) {
		//todo
		// calculate the number of edges each thread will iterate over
		// calculate the number of edges in the wings, tail
		// calculate the number of edges in the dense upper right quadrant (wing width)
		// compute_stats
	}

	uint32_t stripe_len = hyperceiling(wing_width) / 2;
	// evenly divides since both num, denom guaranteed to be powers of 2
	uint32_t n_qs_per_stripe = stripe_len / q_side_len;
	uint32_t right_wing_width = n - wing_width;
	uint32_t n_stripes_in_right_wing = (right_wing_width + stripe_len - 1) / stripe_len;
	fmt::print("{} {} {} {} {}\n", n, m, n_lw_qs, n_rw_qs, n_tail_qs);

	print_arr<uint32_t>(cumulative_n_qs_per_right_wing_stripe, n_stripes_in_right_wing + 1);

	// set the number of threads == number of available cores
//	setenv("OMP_PLACES", "threads", 1);
	setenv("OMP_PLACES", "cores", 1);
//	setenv("OMP_PROC_BIND", "spread", 1);
	setenv("OMP_PROC_BIND", "close", 1);
	int n_cores = omp_get_num_places();
	int partition_n_places = omp_get_partition_num_places();
	int n_threads = n_cores;
	omp_set_num_threads(n_threads);
	fmt::print("n_cores: {}\n", n_cores);
	fmt::print("partition_n_places: {}\n", partition_n_places);
	fmt::print("n_threads: {}\n", n_threads);


	// read out degree file
	std::vector<uint32_t> deg(n);
	std::vector<float> mapped_degs(n);

	read_text_degree_file(out_deg_path, deg);
#pragma omp parallel for
	for (uint32_t i = 0; i < n; ++i) {
		mapped_degs[iso_map[i]] = float(deg[i]);
	}

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


#pragma omp parallel for
	for (uint32_t i = 0; i < n; i++) {
		outgoing_contrib[i] = init_score / mapped_degs[i];
	}

	auto start_time = std::chrono::high_resolution_clock::now();



	for (int iter = 0; iter < num_iterations; iter++) {
//		fmt::print("Iteration: {}\n", iter);
		// left wing
		iterate_left_wing(lw_qs, n_lw_qs, q_side_len, incoming_totals, outgoing_contrib, debug);
		merge_incoming_totals(incoming_totals, incoming_total, 0, wing_width, n_threads);

		// right wing
		iterate_right_wing(rw_qs, n_rw_qs,
		                   wing_width, q_side_len,
		                   incoming_totals, outgoing_contrib,
		                   cumulative_n_qs_per_right_wing_stripe,
		                   n_stripes_in_right_wing,
		                   n_qs_per_stripe,
		                   incoming_total, n, n_threads);
		// tail
		iterate_tail(tail_qs, n_tail_qs,
		             q_side_len,
		             incoming_total, outgoing_contrib, n);
// compute the outgoing contributions for the next iteration
//		print_arr_slice(incoming_total, wing_width, wing_width + 10);

// todo each thread use simd parallelism
//#pragma omp parallel for simd schedule(static)
#pragma omp parallel for schedule(static)
		for (uint32_t i = 0; i < n; i++) {
			// read each threads incoming_total and sum
			scores[i] = base_score + alpha * incoming_total[i];
			outgoing_contrib[i] = scores[i] / mapped_degs[i];
			incoming_total[i] = 0;
		}

	}

	auto end_time = std::chrono::high_resolution_clock::now();

	uint64_t runtime = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	fmt::print("runtime: {}\n", runtime);


#pragma omp parallel for
	for (uint32_t i = 0; i < n; ++i) {
		mapped_results[i] = scores[iso_map[i]];
	}
	bool valid_pr = std::equal(dpl::execution::par_unseq,
	                           mapped_results.begin(), mapped_results.end(),
	                           pr.begin(),
	                           [](double value1, double value2) {
		                           constexpr double epsilon = 1e-4;
		                           return std::fabs(value1 - value2) < epsilon;
	                           });

	uint32_t first_n = 32;
	for (uint32_t i = 0; i < first_n; ++i){
		fmt::print("{} {}\n", mapped_results[i], pr[i]);
	}
	fmt::print("valid_pr: {}\n", valid_pr);

	// clean up
	delete[]lw_qs;
	delete[]rw_qs;
	delete[]tail_qs;
}