//
// Created by atrostan on 13/09/22.
//

#include "PageRank.h"


PageRank::PageRank(int n, int i, std::vector<Edge> &edges, float a) : edges(edges) {
	num_nodes = n;
	num_iters = i;
	alpha = a;
	init();
}

void PageRank::calc_out_degrees() {
	for (auto &e: edges) {
		size_t src = e.source;
		deg[src] += 1.0;
	}
}


void PageRank::init() {
//	src.resize(num_nodes);
//	src.fill(0.0);
//	dst.resize(num_nodes, 0.0);
	deg.resize(num_nodes);
	deg.fill(0.0);
	scores.resize(num_nodes);
	scores.resize(1.0 / num_nodes);
	calc_out_degrees();
}

void PageRank::compute() {
	// measure the total time to complete all PR iterations
	const double init_score = 1.0 / num_nodes;
	const double base_score = (1.0 - alpha) / num_nodes;

	pvector<double> outgoing_contrib(num_nodes);
	pvector<double> incoming_total(num_nodes, 0);

#pragma omp parallel for
	for (uint32_t n = 0; n < num_nodes; n++) {
//		dst[n] = init_score / deg[n];
		outgoing_contrib[n] = init_score / deg[n];
	}
	auto start_time = std::chrono::high_resolution_clock::now();

	for (int iter = 0; iter < num_iters; iter++) {
//		for (int n = 0; n < num_nodes; n++) {
//			src[n] = alpha * dst[n] / deg[n];
//			dst[n] = base_score;
//			incoming_total[n] = 0;
//		}
		for (const auto &e: edges) {
			ul src = e.source;
			ul dst = e.dest;
//			dst[y] += src[x];
			incoming_total[dst] += outgoing_contrib[src];
		}
		for (int n = 0; n < num_nodes; n++) {
			scores[n] = base_score + alpha * incoming_total[n];
			outgoing_contrib[n] = scores[n] / deg[n];
			incoming_total[n] = 0;
		}
	}
	auto end_time = std::chrono::high_resolution_clock::now();
	runtime = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

}


void PageRank::write(std::string path) {
	std::ofstream output_file(path);

	for (int n = 0; n < num_nodes; ++n) {
//		fmt::print("dst[n]: {}\n", dst[n]);
		output_file << std::fixed << std::showpoint;
		output_file << std::setprecision(15);
		output_file << dst[n] << std::endl;
	}

	output_file.close();
}