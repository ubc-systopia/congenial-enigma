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
#include <algorithm>
#include "furhilbert.h"

class QuadrantCalculator {
public:
	int critical_depth = 0;
	int q_idx = 0;
	uint32_t n;
	std::vector<Quadrant> &qs;
	std::vector<std::vector<Quadrant *>> quads;

	QuadrantCalculator(int critical_depth, uint32_t len_n, std::vector<Quadrant> &quads);

	void compute_critical_depth();

	void compute_quadrants();

	void compute_section(int depth, int rot, uint32_t start_x, uint32_t start_y, uint32_t end_x,
	                     uint32_t end_y);


	void compute_directions();

	int get_rot(Corner start, Corner end);



	Direction get_expected_direction(bool flip, Direction hiloop_d);
	Corner get_start_corner(Corner prev_end_corner, Direction prev_dir);
	Corner get_end_corner(Direction d, Corner start);
};


#endif //GRAPH_PREPROCESS_QUADRANTCALCULATOR_H
