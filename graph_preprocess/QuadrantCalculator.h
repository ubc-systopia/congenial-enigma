//
// Created by atrostan on 14/09/22.
//

#ifndef GRAPH_PREPROCESS_QUADRANTCALCULATOR_H
#define GRAPH_PREPROCESS_QUADRANTCALCULATOR_H


#include <vector>
#include "typedefs.h"
#include <math.h>       /* pow */
#include "fmt/core.h"
#include "util.h"

class QuadrantCalculator {
public:
	int n_threads;
	int critical_depth = 0;
	int q_idx = 0;
	uint32_t n; // side length of adjacency matrix
	std::vector<Quadrant> &qs;

	QuadrantCalculator(int critical_depth, uint32_t n, std::vector<Quadrant> &quads);

	void compute_critical_depth();

	void compute_quadrants();

	void compute_section(int depth, int rot, uint32_t start_x, uint32_t start_y, uint32_t end_x,
	                     uint32_t end_y);
};


#endif //GRAPH_PREPROCESS_QUADRANTCALCULATOR_H
