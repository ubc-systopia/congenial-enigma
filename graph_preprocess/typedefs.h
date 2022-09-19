//
// Created by atrostan on 30/08/22.
//

#include <chrono>
#include <string>
#include <sstream>
#ifndef GRAPH_PREPROCESS_TYPEDEFS_H
#define GRAPH_PREPROCESS_TYPEDEFS_H

typedef uint64_t ull;

typedef uint32_t ul;

using namespace std::chrono;
typedef std::chrono::duration<long long, std::milli> time_unit;

enum Order {
	Row,
	Column,
	Hilbert,
//	fgf, // TODO
	End
};

struct PRExptRow {
	std::string graph_name;
	std::string datetime;
	int expt_num;
	int num_iters;
	std::string vertex_order;
	std::string edge_order;
	uint64_t runtime;
};

enum Corner {
	top_right,
	top_left,
	bot_right,
	bot_left,
};

enum Direction {
	left,
	up,
	right,
	down,
};
//
//enum TraversalDirection {
//	upper_right,
//	upper_left,
//	right_down,
//	right_up,
//	bottom_right,
//	bottom_left,
//	left_down,
//	left_up
//};

struct Quadrant {
	uint32_t start_x;
	uint32_t end_x;
	uint32_t start_y;
	uint32_t end_y;
	int rot;
	int idx;
	int hidx;
	Direction d;
	Direction expected_dir;
	Direction intended_dir;
	Corner start;
	Corner end;
};

struct Edge {
	ul source;
	ul dest;
	ull idx;
};

#endif //GRAPH_PREPROCESS_TYPEDEFS_H
