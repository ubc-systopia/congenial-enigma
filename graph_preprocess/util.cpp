//
// Created by atrostan on 30/08/22.
//

#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include "util.h"
#include "omp.h"
#include "fmt/core.h"
#include "fmt/ranges.h"

bool sortByDescendingDegree(const vertex &lhs, const vertex &rhs) {
	if (lhs.degree != rhs.degree)
		return lhs.degree > rhs.degree;
	return lhs.id < rhs.id;
}

void igraph_place_hubs(igraph_t &g, const int k, ul &hub_idx, std::vector<ul> &rank) {
	igraph_vector_t igraph_deg;
	igraph_vector_init(&igraph_deg, 0);
	ul i;

	// select k highest degree vertices
	igraph_degree(&g, &igraph_deg, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
	// use basic hub ordering for now TODO
	std::vector<vertex> vertices(g.n);
	for (ul j = 0; j < g.n; ++j) {
		vertex v;
		v.id = j;
		v.degree = VECTOR(igraph_deg)[j];
		vertices[j] = v;
	}

	std::sort(dpl::execution::par_unseq, vertices.begin(), vertices.end(), sortByDescendingDegree);

	igraph_vector_t hubs_to_remove;
	igraph_vector_init(&hubs_to_remove, k);
	ul j;
	for (j = 0; j < k; ++j) {
		vertex v = vertices[j];
		VECTOR(hubs_to_remove)[j] = v.id;
		rank[hub_idx + j] = VAN(&g, "orig_id", v.id);
	}
	hub_idx += k;

	igraph_delete_vertices(&g, igraph_vss_vector(&hubs_to_remove));

	igraph_vector_destroy(&igraph_deg);
	igraph_vector_destroy(&hubs_to_remove);
}


void create_isomorphism_map(std::vector<ul> &rank, std::map<ul, ul> &map) {
	for (ul i = 0; i < rank.size(); ++i) {
		map[rank[i]] = i;
	}
}

void par_translate_edge_list(std::vector<Edge> &indexed_edges,
                             std::vector<Edge> &mapped_edges,
                             std::vector<ul> &iso_map, ull m) {

//
//	int n_threads = omp_get_max_threads();
//	int tid;
//
//	std::vector<ul> start(n_threads);
//	std::vector<ul> end(n_threads);
//
//	ul n_edges_per_thread = m / n_threads;
//
//	for (int i = 0; i < n_threads; ++i) {
//		start[i] = i * n_edges_per_thread;
//		end[i] = (i + 1) * n_edges_per_thread;
//	}
//
//	// truncate the last section since it may be out of bounds
//	end[n_threads - 1] = m - 1;
//	fmt::print("start: {}\n", start);
//	fmt::print("end: {}\n", end);
//
//	fmt::print("n_threads: {}\n", n_threads);
//	omp_set_num_threads(n_threads);
//#pragma omp parallel private(tid) shared(start, end) default(none)
//	{
//		int tid = omp_get_thread_num();
//		fmt::print("tid: {}\n", tid);
//	}
#pragma omp parallel for default(none) shared(m, mapped_edges, iso_map, indexed_edges) //todo
	for (ull i = 0; i < m; ++i) {
		mapped_edges[i].source = iso_map[indexed_edges[i].source];
		mapped_edges[i].dest = iso_map[indexed_edges[i].dest];
	}

}

void assign_hilbert_keys(std::vector<Edge> &edges, ul n) {
	ull hceil = hyperceiling(n);
#pragma omp parallel for default(none) shared(edges, hceil)// todo
	for (ull i = 0; i < edges.size(); ++i) {
		edges[i].idx = xy2d(hceil, edges[i].source, edges[i].dest);
	}
}


void par_sort_edges(std::vector<Edge> &edges, Order ord, ul n) {
	switch (ord) {
		default:
			break;

		case Column:
			std::sort(dpl::execution::par_unseq, edges.begin(), edges.end(),
			          [](const Edge &lhs, const Edge &rhs) {
				          if (lhs.dest == rhs.dest) {
					          return lhs.source < rhs.source;
				          } else {
					          return lhs.dest < rhs.dest;
				          }
			          });
			break;

		case Row:
			std::sort(dpl::execution::par_unseq, edges.begin(), edges.end(),
			          [](const Edge &lhs, const Edge &rhs) {
				          if (lhs.source == rhs.source) {
					          return lhs.dest < rhs.dest;
				          } else {
					          return lhs.source < rhs.source;
				          }
			          });
			break;

		case Hilbert:
			assign_hilbert_keys(edges, n);
			std::sort(dpl::execution::par_unseq, edges.begin(), edges.end(),
			          [](const Edge &lhs, const Edge &rhs) {
				          return lhs.idx < rhs.idx;
			          });

			break;

//		case Fgf: TODO
//			break;

		case End:
			break;
	}
}

uint32_t next_largest_multiple(uint32_t n, uint32_t critical_depth) {
	assert(critical_depth > 0);
	uint32_t multiple = pow(2, critical_depth);
	return ((n + multiple - 1) / multiple) * multiple;
}


std::ostream &operator<<(std::ostream &os, Direction ec) {
	return os << static_cast<int>(ec);
}

void print_quad(Quadrant &q) {
	fmt::print(
			"[idx: {:<4} rot: {:<4} hidx: {:<4} expected: {:<5}] || sx: {:<20} | ex: {:<20} | sy: {:<20} | ey: {:<20} || {:<15}->{:<15}\n",
			q.idx, q.rot, q.hidx, dir_str(q.expected_dir), q.start_x, q.end_x, q.start_y, q.end_y, corner_str(q.start),
			corner_str(q.end));
}

void print_seperator() {
	std::ostringstream os;
	for (int i = 0; i < 86; ++i) {
		os << "#";
	}
	os << "\n";
	fmt::print("{}\n", os.str());
}

std::string dir_str(Direction d) {
	switch (d) {
		case right:
			return "Right";
		case down:
			return "Down";
		case left:
			return "Left";
		case up:
			return "Up";
	}
}

std::pair<uint32_t, uint32_t> rotate_point(uint32_t cx, uint32_t cy, int angle, uint32_t x, uint32_t y) {

	int s;
	int c;
	uint32_t xnew;
	uint32_t ynew;
	if (angle == 0) return std::make_pair(x, y);
	else {
		if (angle == 90) {
			c = 0;
			s = 1;
			xnew = ((x - cx) * c) - ((cy - y) * s) + cx;
			ynew = cy - ((cy - y) * c) + ((x - cx) * s);
		} else if (angle == 180) {
			c = -1;
			s = 0;

			xnew = ((x - cx) * c) - ((cy - y) * s) + cx - 1;
			ynew = cy - ((cy - y) * c) + ((x - cx) * s) - 1;
		} else if (angle == 270) {
			c = 0;
			s = -1;
			xnew = ((x - cx) * c) - ((cy - y) * s) + cx - 1;
			ynew = cy - ((cy - y) * c) + ((x - cx) * s) - 1;
		}
	}


	// translate point back to origin:
//	x -= cx;
//	y -= cy;

	// rotate point

	// translate point back:
//	x = xnew + cx;
//	y = cy - ynew;
	return std::make_pair(xnew, ynew);
}

std::string corner_str(Corner c) {
	switch (c) {
		case top_right:
			return "Top-Right";
		case top_left:
			return "Top-Left";
		case bot_right:
			return "Bottom-Right";
		case bot_left:
			return "Bottom-Left";
	}
}


uint32_t int_log(int base, uint32_t x) {
	return (int) (log(x) / log(base));
}