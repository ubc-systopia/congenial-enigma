//
// Created by atrostan on 13/09/22.
//

#ifndef GRAPH_PREPROCESS_PAGERANK_H
#define GRAPH_PREPROCESS_PAGERANK_H


#include "fmt/core.h"
#include "typedefs.h"
#include <fmt/ranges.h>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <vector>
#include "pvector.h"

class PageRank {
public:
	int num_iters;
	int num_nodes;
	float alpha;
	uint64_t runtime;
	pvector<double> src;
	pvector<double> dst;
	pvector<double> deg;
	pvector<double> scores;


	std::vector<Edge> &edges;

	PageRank(int n, int i, std::vector<Edge> &edges, float alpha);

	void init();

	void calc_out_degrees();

	void compute();

	void write(std::string path);

};


#endif //GRAPH_PREPROCESS_PAGERANK_H
