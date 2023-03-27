//
// Created by atrostan on 30/08/22.
//

#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <eigen3/Eigen/SparseCore>
#include "util.h"
#include "omp.h"
#include "fmt/core.h"
#include "fmt/ranges.h"
#include "pvector.h"


uint32_t get_min_id_in_range(uint32_t u, uint32_t v_min, uint32_t v_max, CSR &out_csr) {
	uint32_t start = out_csr.index[u];
	uint32_t end = out_csr.index[u + 1];
	if (start == end) { return -1; }

	uint32_t min_neigh = out_csr.neighbours[start];
	uint32_t max_neigh = out_csr.neighbours[end - 1];
	uint64_t offset = lower_bound<uint32_t, uint64_t>(
		v_min,
		out_csr.neighbours,
		start,
		end
	);
	uint32_t candidate = out_csr.neighbours[offset];
	if (offset == end) {
		if (max_neigh < v_min) {
			return -1;
		}
	}
	return candidate;
}

/**
 * Get the number of out-neighbours of u, whose ID lies between
 * [v_min, v_max)
 * @param u
 * @param v_min
 * @param v_max
 * @param out_csr - the out-csr representing the graph storing u and its out-neighbours
 * @return
 */
uint32_t n_out_neighbours_between(uint32_t u, uint32_t v_min, uint32_t v_max, CSR &out_csr) {
	uint32_t start = out_csr.index[u];
	uint32_t end = out_csr.index[u + 1];
	uint32_t min_neigh = out_csr.neighbours[start];
	uint32_t max_neigh = out_csr.neighbours[end - 1];
//	uint64_t min_offset, max_offset;
//	auto p1 = binary_search(out_csr.neighbours, v_min, start, end);
	uint64_t min_offset = lower_bound<uint32_t, uint64_t>(
		v_min,
		out_csr.neighbours,
		start,
		end
	);
	uint64_t max_offset = lower_bound<uint32_t, uint64_t>(
		v_max,
		out_csr.neighbours,
		start,
		end
	);
	return max_offset - min_offset;
}

uint32_t filter_empty_quads(Quad *&qs, uint32_t n_quads) {
	uint32_t n_empty = 0;
#pragma omp parallel for reduction(+:n_empty)
	for (uint32_t i = 0; i < n_quads; ++i) {
		if (qs[i].nnz == 0) {
			n_empty++;
		}
	}

	uint32_t reduced_size = n_quads - n_empty;
	Quad *new_qs = new Quad[reduced_size];
	// create a new quad array, copy over the relevant data, and swap the pointers
	uint32_t non_empty_idx = 0;
	for (uint32_t i = 0; i < n_quads; ++i) {
		if (qs[i].nnz == 0) { continue; }
		new_qs[non_empty_idx].qx = qs[i].qx;
		new_qs[non_empty_idx].qy = qs[i].qy;
		new_qs[non_empty_idx].q_idx = qs[i].q_idx;
		new_qs[non_empty_idx].nnz = qs[i].nnz;
		new_qs[non_empty_idx].edges = qs[i].edges;

		qs[i].qx = -1;
		qs[i].qy = -1;
		qs[i].q_idx = -1;
		qs[i].nnz = -1;
		qs[i].edges = nullptr;

		non_empty_idx++;
	}
	qs = new_qs;
	return reduced_size;
}

/**
 *
 * Reorder the edges contained in a quadrant using the hilbert order
	* @param q
	* @param side_len
 */
void horder_edges_in_q(Quad &q, uint32_t side_len) {
	// create a temporary vector to store the reordered edges
	uint32_t nnz = q.nnz;
	struct IndexedEdge {
		uint32_t src;
		uint32_t dest;
		uint32_t h_idx;
	};
	std::vector<IndexedEdge> copy(nnz);
	for (uint32_t j = 0; j < nnz; ++j) {
		uint32_t src = q.edges[j * 2];
		uint32_t dest = q.edges[j * 2 + 1];
		copy[j].src = src;
		copy[j].dest = dest;
		copy[j].h_idx = xy2d(side_len, src, dest);
	}
	std::sort(
		copy.begin(),
		copy.end(),
		[](const IndexedEdge &a, const IndexedEdge &b) -> bool { return a.h_idx < b.h_idx; }
	);
	// copy the (now sorted) edges back to the edges array
	for (uint32_t j = 0; j < nnz; ++j) {
		q.edges[j * 2] = copy[j].src;
		q.edges[j * 2 + 1] = copy[j].dest;
	}
}

/**
 * Given an edge (u, v) that resides in a quadrant whose quadrant side length is q_side_len compute the local (offset) version of that edge
 * @param u
 * @param v
 * @param q_side_len
 * @return
 */
std::pair<uint32_t, uint32_t> edge_offset(uint32_t u, uint32_t v, uint32_t q_side_len) {
	return {u % q_side_len, v % q_side_len};
}

bool sortByDescendingDegree(const vertex &lhs, const vertex &rhs) {
	if (lhs.degree != rhs.degree)
		return lhs.degree > rhs.degree;
	return lhs.id < rhs.id;
}

void igraph_place_hubs(igraph_t &g, const int k, ul &hub_idx, std::vector<ul> &rank) {
	igraph_vector_int_t igraph_deg;
	igraph_vector_int_init(&igraph_deg, 0);
	ul i;

	// select k highest degree vertices
	igraph_degree(&g, &igraph_deg, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
	// use basic hub ordering for now TODO
	std::vector<vertex> vertices(g.n);
	for (ul j = 0; j < g.n; ++j) {
		vertex v;
		v.id = j;

		// igraph counts duplicate neighbours for directed graphs that were
		// ingested as undirected
		// e.g. [0 -> 1, 1 -> 0] --> deg(0) = 2
		// instead iterate over in + outneighbourhoods, and count all the UNIQUE
		// neighbours to compute a vertex's degree
		v.degree = VECTOR(igraph_deg)[j];
		vertices[j] = v;
	}


	std::sort(dpl::execution::par_unseq, vertices.begin(), vertices.end(), sortByDescendingDegree);

//	fmt::print("vertices: \n");
//	for (const auto &v: vertices) {
//		fmt::print("[{} {}], ", v.id, v.degree);
//	}
//	fmt::print("\n");
	igraph_vector_int_t hubs_to_remove;
	igraph_vector_int_init(&hubs_to_remove, k);
	ul j;
	for (j = 0; j < k; ++j) {
		vertex v = vertices[j];
		VECTOR(hubs_to_remove)[j] = v.id;
		rank[hub_idx + j] = VAN(&g, "orig_id", v.id);
	}
	hub_idx += k;

	igraph_delete_vertices(&g, igraph_vss_vector(&hubs_to_remove));

	igraph_vector_int_destroy(&igraph_deg);
	igraph_vector_int_destroy(&hubs_to_remove);
}


void create_isomorphism_map(std::vector<ul> &rank, std::map<ul, ul> &map) {
	for (ul i = 0; i < rank.size(); ++i) {
		map[rank[i]] = i;
	}
}

/**
 * Normalize a PageRank vector so that it can be compared against ground truth
 * @param orig
 * @param normalized
 */
void normalize(std::vector<double> &orig, std::vector<double> &normalized) {
	double sum = 0.0;
//#pragma omp parallel for reduction(+:sum)
	for (int n = 0; n < orig.size(); n++) {
		sum += orig[n];
	}

//#pragma omp parallel for
	for (int n = 0; n < orig.size(); n++) {
		normalized[n] = orig[n] / sum;
	}
}

void par_translate_edge_list(std::vector<Edge> &indexed_edges,
                             std::vector<Edge> &mapped_edges,
                             std::vector<ul> &iso_map, ull m) {

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

bool is_power_of_2(uint32_t n) {
	return (n > 0 && ((n & (n - 1)) == 0));
}
