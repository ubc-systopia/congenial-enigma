//
// Created by atrostan on 14/09/22.
//

#include <cassert>
#include "QuadrantCalculator.h"


bool compare_quad_by_x_y(const Quadrant &a, const Quadrant &b) {
	return a.start_y < b.start_y ||
	       a.start_y == b.start_y && a.start_x < b.start_x;
}


bool compare_quad_by_hidx(const Quadrant &a, const Quadrant &b) {
	return a.hidx < b.hidx;
}

Corner QuadrantCalculator::get_start_corner(Corner prev_end_corner, Direction prev_dir) {
	if (prev_end_corner == top_right) {
		if (prev_dir == right) return top_left;
		if (prev_dir == up) return bot_right;
//		if (prev_dir == down) fmt::print("henlo\n");
	}
	if (prev_end_corner == bot_left) {
		if (prev_dir == down) return top_left;
		if (prev_dir == left) return bot_right;
	}
}

Corner QuadrantCalculator::get_end_corner(Direction d, Corner start) {
	if (start == top_left) {
		if (d == right || d == up) return top_right;
		if (d == down || d == left) return bot_left;
	}
	if (start == bot_right) {
		if (d == up || d == right) return top_right;
		if (d == down || d == left) return bot_left;
	} else {
		fmt::print("ERROR: {} {} \n", dir_str(d), corner_str(start));
	}
}

int QuadrantCalculator::get_rot(Corner start, Corner end) {
	if      (start == top_left && end == top_right) return 0;
	else if (start == top_left && end == bot_left)  return 90;
	else if (start == bot_right && end == bot_left) return 180;
	else if (start == bot_right && end == top_right) return 270; // todo reversed
	else {
		return 0;
	}
}

QuadrantCalculator::QuadrantCalculator(int cd, uint32_t len_n, std::vector<Quadrant> &qds) : qs(qds) {
	n = len_n;
	critical_depth = cd;
	compute_quadrants();
	// sort quadrants by row and column
	std::sort(qs.begin(), qs.end(), compare_quad_by_x_y);
	int qlen = int(sqrt(qs.size()));

	// the directions given by HILLOOP_d need to be flipped when the
	// side length of the quadrant matrix is an ODD power of four
	bool flip_hiloop_d = int_log(4, qs.size()) % 2 == 1;

	fmt::print("flip_hiloop_d: {}\n", flip_hiloop_d);
	
	quads.resize(qlen);
	for (int i = 0; i < qlen; ++i) {
		quads[i].resize(qlen);
	}

	int k = 0;
	for (int i = 0; i < qlen; ++i) {
		for (int j = 0; j < qlen; ++j) {
			quads[i][j] = &qs[k];
			++k;
		}
	}
	// assign each quadrant in the quadrant matrix its furhilbert index
	int qi, qj;
	qi = qj = k = 0;
	struct Pt {
		int x;
		int y;
	};
	int prev_qi = 0;
	int prev_qj = 0;
	FUR_HILBERT_FOR(qi, qj, 0, qlen, 0, qlen)
						{

							int HILLOOP_d = HILLOOP_nanoL & 1 | HILLOOP_nanoH & 2;
							Direction expected = get_expected_direction(flip_hiloop_d, Direction(HILLOOP_d));
							if (flip_hiloop_d) {
								quads[qi][qj]->expected_dir = expected;
								quads[qi][qj]->hidx = k;
							} else {
								quads[qj][qi]->expected_dir = expected;
								quads[qj][qi]->hidx = k;
							}

							++k;
						}
	FUR_HILBERT_END(qi, qj);
	std::sort(qs.begin(), qs.end(), compare_quad_by_hidx);
	qs[0].start = top_left;
	qs[0].end = top_right;
	Corner prev_end_corner = top_right;
	Direction prev_dir = qs[0].expected_dir;
	for (int i = 1; i < qs.size(); ++i) {
		qs[i].start = get_start_corner(prev_end_corner, prev_dir);
		qs[i].end = get_end_corner(qs[i].expected_dir, qs[i].start);
		prev_end_corner = qs[i].end;
		prev_dir = qs[i].expected_dir;
	}
	for (int i = 0; i < qs.size(); ++i) {
		qs[i].rot = get_rot(qs[i].start, qs[i].end);
	}
}

void QuadrantCalculator::compute_quadrants() {
	compute_section(0, 0, 0, 0, n, n);
}

Direction QuadrantCalculator::get_expected_direction(bool flip, Direction hiloop_d) {
	if (flip) {
		switch (hiloop_d) {
			case right:
				return down;
			case down:
				return right;
			case up:
				return left;
			case left:
				return up;
		}
	} else {
		return hiloop_d;
	}
}

void QuadrantCalculator::compute_directions() {
	// compute the seed direction
	Quadrant &seed = qs[0];
//	seed.expected_dir = get_expected_direction(seed.end_x - seed.start_x, seed.end_y - seed.start_y);
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

//		assert(end_x - start_x == end_y - start_y);

		qs[q_idx] = Quadrant{start_x, end_x, start_y, end_y, rot, q_idx, 0, right};
		++q_idx;
	} else {
		compute_section(depth + 1, rot, start_x, start_y, mid_x, mid_y); // upper left
		compute_section(depth + 1, rot, mid_x, start_y, end_x, mid_y); // upper right
		compute_section(depth + 1, rot, start_x, mid_y, mid_x, end_y); // lower left
		compute_section(depth + 1, rot, mid_x, mid_y, end_x, end_y); // lower right
	}
}
