//
// Created by atrostan on 30/08/22.
//

#include <chrono>
#include <string>

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

struct Quadrant {
	uint32_t start_x;
	uint32_t end_x;
	uint32_t start_y;
	uint32_t end_y;
	int rot;
	int idx;
};

struct Edge {
	ul source;
	ul dest;
	ull idx;
};

#endif //GRAPH_PREPROCESS_TYPEDEFS_H
