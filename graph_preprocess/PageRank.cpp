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
	src.resize(num_nodes, 0.0);
	dst.resize(num_nodes, 0.0);
	deg.resize(num_nodes, 0.0);

	calc_out_degrees();
}

void PageRank::compute() {
	// measure the total time to complete all PR iterations
	auto start_time = std::chrono::high_resolution_clock::now();

	for (int iter = 0; iter < num_iters; iter++) {
		for (int n = 0; n < num_nodes; n++) {
			src[n] = alpha * dst[n] / deg[n];
			dst[n] = 1.0 - alpha;
		}
		for (const auto &e: edges) {
			ul x = e.source;
			ul y = e.dest;

			dst[y] += src[x];
		}
	}
	auto end_time = std::chrono::high_resolution_clock::now();
	runtime = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

	// normalize
	float sum = 0.0;
	for (int n = 0; n < num_nodes; n++) {
		sum += dst[n];
	}

	for (int n = 0; n < num_nodes; n++) {
		dst[n] = dst[n] / sum;
	}

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