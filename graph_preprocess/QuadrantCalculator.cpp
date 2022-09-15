//
// Created by atrostan on 14/09/22.
//

#include "QuadrantCalculator.h"

QuadrantCalculator::QuadrantCalculator(int cd, uint32_t side_len, std::vector<Quadrant> &quads) : qs(quads) {
	n = side_len;
	critical_depth = cd;
	compute_quadrants();
}

void QuadrantCalculator::compute_quadrants() {
	compute_section(0, 0, 0, 0, n - 1, n - 1);
}


void QuadrantCalculator::compute_section(int depth, int rot, uint32_t start_x, uint32_t start_y, uint32_t end_x,
                                         uint32_t end_y) {
	uint32_t len_x = end_x - start_x;
	uint32_t len_y = end_y - start_y;

	uint32_t dec_x = (len_x) % 2;
	uint32_t dec_y = (len_y) % 2;

//	uint32_t mid_x = ceil((start_x * 1.0 + end_x) / 2) - dec_x;
//	uint32_t mid_y = ceil((start_y * 1.0 + end_y) / 2) - dec_y;

	uint32_t mid_x = (start_x + end_x) / 2;
	uint32_t mid_y = (start_y + end_y) / 2;

	if (depth == critical_depth) {
		qs[q_idx] = Quadrant{start_x, end_x, start_y, end_y, rot, q_idx};
		++q_idx;
	} else {
		compute_section(depth + 1, rot, start_x, start_y, mid_x - 1, mid_y - 1); // upper left
		compute_section(depth + 1, rot, mid_x, start_y, end_x, mid_y - 1); // upper right
		compute_section(depth + 1, rot, start_x, mid_y, mid_x - 1, end_y); // lower left
		compute_section(depth + 1, rot, mid_x, mid_y, end_x, end_y); // lower right

	}
}
