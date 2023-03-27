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

    uint64_t size() const { return _n_threads * _n_vertices * sizeof(float); }

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

            incoming_total[offset + u] += grid(t, u);
            grid(t, u) = 0;
        }


}
//
//void right_wing_tasks(Quad *qs,
//                      uint32_t wing_width, uint32_t q_side_len,
//                      float *outgoing_contrib,
//                      uint32_t *cumulative_n_qs_per_right_wing_stripe,
//                      uint32_t n_stripes_in_right_wing,
//                      uint32_t n_qs_per_stripe,
//                      float *incoming_total,
//                      uint32_t n, int n_cores) {
//    for (uint32_t i = 0; i < n_stripes_in_right_wing; ++i) {
//        uint32_t stripe_start = cumulative_n_qs_per_right_wing_stripe[i];
//        uint32_t stripe_end = cumulative_n_qs_per_right_wing_stripe[i + 1];
//        uint32_t v_start = wing_width + (n_qs_per_stripe * q_side_len * i);
//        uint32_t v_end = wing_width + (n_qs_per_stripe * q_side_len * (i + 1));
//        v_end = (v_end > n) ? n : v_end;
//        uint32_t v_len = v_end - v_start;
////https://stackoverflow.com/questions/13065943/task-based-programming-pragma-omp-task-versus-pragma-omp-parallel-for
//        // spawn tasks for right wing
//#pragma omp taskgroup task_reduction(+:incoming_total[v_start:v_len])
//        {
//            for (uint32_t j = stripe_start; j < stripe_end; ++j) {
//                uint32_t task_priority = stripe_end - j;
//                Quad &q = qs[j];
//                uint32_t qx = q_side_len * q.qx;
//                uint32_t qy = q_side_len * q.qy;
//#pragma omp task shared(outgoing_contrib) in_reduction(+:incoming_total[v_start:v_len])  priority(task_priority)
//                {
//                    Quad &q = qs[j];
//                    for (uint32_t k = 0; k < q.nnz; ++k) {
//                        uint32_t src = q.edges[k * 2];
//                        uint32_t dest = q.edges[k * 2 + 1];
//                        uint32_t u = q_side_len * q.qx + src;
//                        uint32_t v = q_side_len * q.qy + dest + wing_width;
//                        incoming_total[v] += outgoing_contrib[u];
//                    }
//                }
//            }
//        }
//    }
//}



void right_wing_tasks(Quad *qs,
                      uint32_t wing_width, uint32_t q_side_len,
                      float *outgoing_contrib,
                      uint32_t *cumulative_n_qs_per_right_wing_stripe,
                      uint32_t n_stripes_in_right_wing,
                      uint32_t n_qs_per_stripe,
                      float *incoming_total,
                      uint32_t n, int n_cores) {
	for (uint32_t i = 0; i < n_stripes_in_right_wing; ++i) {
		uint32_t stripe_start = cumulative_n_qs_per_right_wing_stripe[i];
		uint32_t stripe_end = cumulative_n_qs_per_right_wing_stripe[i + 1];
		uint32_t v_start = wing_width + (n_qs_per_stripe * q_side_len * i);
		uint32_t v_end = wing_width + (n_qs_per_stripe * q_side_len * (i + 1));
		v_end = (v_end > n) ? n : v_end;
		uint32_t v_len = v_end - v_start;
//https://stackoverflow.com/questions/13065943/task-based-programming-pragma-omp-task-versus-pragma-omp-parallel-for
		// spawn tasks for right wing
#pragma omp parallel num_threads(n_cores)
#pragma omp single nowait
		{
#pragma omp taskgroup task_reduction(+:incoming_total[v_start:v_len])
			{
				for (uint32_t j = stripe_start; j < stripe_end; ++j) {
					uint32_t task_priority = stripe_end - j;
					Quad &q = qs[j];
					uint32_t qx = q_side_len * q.qx;
					uint32_t qy = q_side_len * q.qy;
#pragma omp task shared(outgoing_contrib) \
        in_reduction(+:incoming_total[v_start:v_len]) \
        priority(task_priority) \
        affinity(outgoing_contrib[qx:qx+q_side_len], incoming_total[qy:qy+q_side_len])
					{
						Quad &q = qs[j];
						for (uint32_t k = 0; k < q.nnz; ++k) {
							uint32_t src = q.edges[k * 2];
							uint32_t dest = q.edges[k * 2 + 1];
							uint32_t u = q_side_len * q.qx + src;
							uint32_t v = q_side_len * q.qy + dest + wing_width;
							incoming_total[v] += outgoing_contrib[u];
						}
					}
				}
			}
		}
	}
}

void iterate_left_wing(Quad *qs, uint32_t n_qs,
                       uint32_t q_side_len,
                       PageRankGrid &grid, float *outgoing_contrib,
                       float *incoming_total,
                       uint32_t wing_width,
                       bool debug,
                       int n_cores) {
//	std::vector<uint16_t> thread_assignments;
//	 if debugging, keep track of which thread computed which quadrant
//	if (debug) {
//		thread_assignments.resize(n_qs);
//	}
#pragma omp parallel num_threads(n_cores)
    {
        int tid = omp_get_thread_num();
//#pragma omp for schedule(runtime)
//#pragma omp for schedule(static, 1)
//#pragma omp for schedule(dynamic, 1)
//#pragma omp for schedule(static, 1) reduction(+:incoming_total[:wing_width])
#pragma omp for schedule(dynamic, 1) reduction(+:incoming_total[:wing_width])
        for (uint32_t i = 0; i < n_qs; ++i) {
            Quad &q = qs[i];
//			if (debug) {
//				thread_assignments[i] = tid;
//			}
            //todo would the following benefit from vectorization? (test)
//#pragma omp simd
            for (uint32_t j = 0; j < q.nnz; ++j) {
                uint32_t src = q.edges[j * 2];
                uint32_t dest = q.edges[j * 2 + 1];
                uint32_t u = q_side_len * q.qx + src;
                uint32_t v = q_side_len * q.qy + dest;
//				grid[tid][v] += outgoing_contrib[u];
//				grid(tid, v) += outgoing_contrib[u];
                incoming_total[v] += outgoing_contrib[u];
            }
        }
    }
//	if (debug) {
//		fmt::print("thread_assignments: {}\n", thread_assignments);
//		omp_sched_t sched;
//		int chunk;
//		omp_get_schedule(&sched, &chunk);
//		fmt::print("sched, chunk: {} {}\n", sched, chunk);
//	}


}

void iterate_tail(Quad *qs, uint32_t n_qs,
                  uint32_t q_side_len,
                  float *incoming_total, float *outgoing_contrib,
                  uint32_t n, int num_cores) {
    // each rectangle in the tail is guaranteed to write to an exclusive
    // range of destination ids
    // rectangles are sorted by qy
//#pragma omp parallel for schedule(static, 1) num_threads(num_cores)
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_cores)
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
                        uint32_t n, uint16_t n_threads, int num_cores) {
#pragma omp parallel num_threads(num_cores)
    {
        int tid = omp_get_thread_num();
        for (uint32_t i = 0; i < n_stripes_in_right_wing; ++i) {
            uint32_t stripe_start = cumulative_n_qs_per_right_wing_stripe[i];
            uint32_t stripe_end = cumulative_n_qs_per_right_wing_stripe[i + 1];
//			uint32_t u_start = wing_width + (n_qs_per_stripe * i);
//			uint32_t u_end = wing_width + (n_qs_per_stripe * (i + 1));
            uint32_t v_start = wing_width + (n_qs_per_stripe * q_side_len * i);
            uint32_t v_end = wing_width + (n_qs_per_stripe * q_side_len * (i + 1));
            v_end = (v_end > n) ? n : v_end;
            uint32_t v_len = v_end - v_start;
//#pragma omp for schedule(static, 1)
//#pragma omp for schedule(static, 1) reduction(+:incoming_total[v_start:v_len])
#pragma omp for schedule(dynamic, 1) reduction(+:incoming_total[v_start:v_len])
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
//					uint32_t v = q_side_len * qy + dest;
                    uint32_t v = q_side_len * q.qy + dest + wing_width;
//					fmt::print("v_start, v_end: {}  {}\n", v_start, v_end);
//					fmt::print("tid, u, v: {} {} {} {} {} \n", tid, u, v, v_start, v_end);

                    incoming_total[v] += outgoing_contrib[u];
//					grid(tid, v) += outgoing_contrib[u];
//					grid[tid][v] += outgoing_contrib[u];
                }
            }

            // once we've completed this stripe, merge the results into
            // the global incoming total array
//			uint32_t v_offset = wing_width + (n_qs_per_stripe * q_side_len * i);
//			uint32_t v_end = wing_width + (n_qs_per_stripe * q_side_len * (i + 1));
//			v_end = (v_end > n) ? n : v_end;
//			v_end -= v_offset;
//			merge_incoming_totals(grid, incoming_total, v_offset, v_end, n_threads);

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
    int n_cores = 0;
    bool debug = false;
    std::string schedule = ""; // 1 of {"static", "dynamic"}

    while ((opt = getopt(argc, argv, "bg:d:l:e:w:i:s:c:")) != -1) {
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
            case 'c':
                n_cores = atoi(optarg);
                break;
        }
    }
    assert(is_power_of_2(q_side_len));
//    fmt::print("_OPENMP: {}\n", _OPENMP);
    fmt::print("KMP_VERSION_MAJOR: {}\n", KMP_VERSION_MAJOR);
    fmt::print("n_cores: {}\n", n_cores);

//	return 0;

    fmt::print("kmp: {}\n", kmp_get_stacksize());
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
    fmt::print("Read left wing.\n");
    n_lw_qs = res[0];
    n = res[1];
    m = res[2];
    res = read_quad_array(rw_path, rw_qs, false, true, cumulative_n_qs_per_right_wing_stripe);
    fmt::print("Read right wing.\n");

    n_rw_qs = res[0];
    res = read_quad_array(tail_path, tail_qs, true, false, empty);
    fmt::print("Read tail wing.\n");

    n_tail_qs = res[0];

    for (uint32_t i = 0; i <10; ++i) {
        fmt::print("left_wing_qs[i].q_idx: {}\n", lw_qs[i].qx);
        fmt::print("right_wing_qs[i].q_idx: {}\n", lw_qs[i].qy);
    }

    if (debug) {
        //todo
        // calculate the number of edges each thread will iterate over
        // calculate the number of edges in the wings, tail
        // calculate the number of edges in the dense upper right quadrant (wing width)
        // compute_stats
    }

    uint32_t stripe_len = hyperceiling(wing_width) / 2;
    fmt::print("stripe_len: {}\n", stripe_len);
    // evenly divides since both num, denom guaranteed to be powers of 2
    uint32_t n_qs_per_stripe = stripe_len / q_side_len;
    uint32_t right_wing_width = n - wing_width;
    uint32_t n_stripes_in_right_wing = (right_wing_width + stripe_len - 1) / stripe_len;
    fmt::print("{} {} {} {} {}\n", n, m, n_lw_qs, n_rw_qs, n_tail_qs);

    print_arr<uint32_t>(cumulative_n_qs_per_right_wing_stripe, n_stripes_in_right_wing + 1);

    // set the number of threads == number of available cores (e.g. cores, threads)
//	setenv("OMP_PLACES", "cores", 1);
//	setenv("OMP_PROC_BIND", "close", 1);

    int default_omp_stack_size = kmp_get_stacksize() ; // the amount of stack memory assigned to each thread, in bytes
    uint32_t required_stack_size = wing_width * 4;

    uint32_t stack_size = hyperceiling(required_stack_size);
    fmt::print(" {} {} {}\n", default_omp_stack_size, required_stack_size, stack_size);
    fmt::print("kmp_get_stacksize() / 1024: {}\n", kmp_get_stacksize());
    if (stack_size > default_omp_stack_size) {
        std::string stack_size_str = fmt::format("{}B", stack_size);
        fmt::print("Setting OMP_STACKSIZE to {}.\n",  stack_size_str);
//		setenv("OMP_STACKSIZE", stack_size_str.c_str(), 1);
        setenv("OMP_STACKSIZE", "25m", 1);
    }
    fmt::print("kmp_get_stacksize() / 1024: {}\n", kmp_get_stacksize());

    int omp_n_cores = omp_get_num_places();
    int partition_n_places = omp_get_partition_num_places();
    int n_threads = omp_n_cores;
    if (n_cores == 0) n_cores = n_threads;
    fmt::print("omp_get_max_threads(): {}\n", omp_get_max_threads());

//    omp_set_num_threads(n_threads);
//    fmt::print("omp_n_cores: {}\n", omp_n_cores);
//    fmt::print("partition_n_places: {}\n", partition_n_places);
//    fmt::print("n_threads: {}\n", n_threads);
//	uint16_t n_threads = omp_get_max_threads();
//	fmt::print("n_threads: {}\n", n_threads);



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
//	PageRankGrid incoming_totals = PageRankGrid(n_threads, wing_width);
    PageRankGrid incoming_totals = PageRankGrid(0, 0);
    fmt::print("incoming_totals.size(): {}\n", incoming_totals.size());

//	return 0;
    // global arrays
    float *incoming_total = new float[n]();
    float *outgoing_contrib = new float[n]();
    float *scores = new float[n]();


#pragma omp parallel for
    for (uint32_t i = 0; i < n; i++) {
        outgoing_contrib[i] = init_score / mapped_degs[i];
    }

    auto start_time = std::chrono::high_resolution_clock::now();


    fmt::print("size_of(float): {}\n", sizeof(float));
    for (int iter = 0; iter < num_iterations; iter++) {
//		fmt::print("Iteration: {}\n", iter);
        // left wing
        iterate_left_wing(lw_qs, n_lw_qs, q_side_len, incoming_totals, outgoing_contrib, incoming_total, wing_width, debug, n_cores);
//		merge_incoming_totals(incoming_totals, incoming_total, 0, wing_width, n_threads);

        // right wing
        iterate_right_wing(rw_qs, n_rw_qs,
                           wing_width, q_side_len,
                           incoming_totals, outgoing_contrib,
                           cumulative_n_qs_per_right_wing_stripe,
                           n_stripes_in_right_wing,
                           n_qs_per_stripe,
                           incoming_total, n, n_threads, n_cores);
        // tail
        iterate_tail(tail_qs, n_tail_qs,
                     q_side_len,
                     incoming_total, outgoing_contrib, n, n_cores);
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