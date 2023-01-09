//
// Created by atrostan on 28/12/22.
//

#include "hilburn_preproc.h"

struct pair_hash {
	template<class T1, class T2>
	size_t operator()(const std::pair<T1, T2> &pair) const {
		std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
	}
};

typedef std::unordered_map<std::pair<uint32_t, uint32_t>, bool, pair_hash> edge_list_map;

bool sort_by_q_idx(Quad *x, Quad *y) { return x->q_idx < y->q_idx; }

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

void print_out_neighs(uint32_t u, CSR &out_csr) {
	if (u > out_csr.num_nodes) return;
	uint32_t start = out_csr.index[u];
	uint32_t end = out_csr.index[u + 1];
	fmt::print("{} --> ", u);
	for (uint32_t offset = start; offset < end; ++offset) {
		fmt::print("{}, ", out_csr.neighbours[offset]);
	}
	fmt::print("\n");
}

template<typename T>
void print_arr(T *arr, uint64_t size) {
	fmt::print("[");
	for (uint64_t i = 0; i < size - 1; ++i) { fmt::print("{},", arr[i]); }
	fmt::print("{}]\n", arr[size - 1]);
}

uint32_t n_out_neighbours(uint32_t u, CSR &out_csr) {
	if (u > out_csr.num_nodes) return 0;
	return out_csr.index[u + 1] - out_csr.index[u];
}

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

uint32_t get_max_id_in_range(uint32_t u, uint32_t v_min, uint32_t v_max, CSR &out_csr) {
	uint32_t start = out_csr.index[u];
	uint32_t end = out_csr.index[u + 1];
	if (start == end) { return 0; }
	uint32_t min_neigh = out_csr.neighbours[start];
	uint32_t max_neigh = out_csr.neighbours[end - 1];
//	uint64_t offset = lower_bound<uint32_t, uint64_t>(
//		v_max,
//		out_csr.neighbours,
//		start,
//		end
//	);
	uint64_t offset = binary_search(out_csr.neighbours, v_max, start, end).second;
	uint32_t candidate = out_csr.neighbours[offset];
//	print_out_neighs(u, out_csr);
//	fmt::print("{} {} {} {}\n", offset, end, candidate, v_max);
	if (offset == end) {
		if (max_neigh <= v_max && max_neigh >= v_min) {
//			fmt::print("max_neigh: {}\n", max_neigh);
			return max_neigh;
		} else {
			return 0;
		}
	} else if (offset < end) {
		if (candidate > v_max) {
			if (out_csr.neighbours[offset - 1] < v_max) return out_csr.neighbours[offset - 1];
			else {
				return 0;
			}
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


void read_degs(std::string out_path, std::string in_path, std::vector<uint32_t> &out_degs,
               std::vector<uint32_t> &in_degs, uint32_t n, std::vector<uint32_t> &iso_map) {
	std::vector<uint32_t> ods(n);
	std::vector<uint32_t> ids(n);

	read_text_degree_file(out_path, ods);
	read_text_degree_file(in_path, ids);

#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < n; ++i) {
		uint32_t parsb_id = iso_map[i];
		out_degs[parsb_id] = ods[i];
		in_degs[parsb_id] = ids[i];
	}

}

uint32_t mod(const uint32_t v, const uint32_t n) {
//	assert(is_power_of_2(n));
	return v & (n - 1);
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
//	return {mod(u, q_side_len), mod(v, q_side_len)};
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


bool check_edge(uint32_t src, uint32_t dest, edge_list_map &map) {
	std::pair<uint32_t, uint32_t> key = {src, dest};
//	fmt::print("key, map[key]: {} {}\n", key, map[key]);
	// incorrect mapping
	if (map.find(key) == map.end()) { return false; }
	else {
		bool in_map = map[key];
		// already seen this edge
		if (in_map) { return false; }
		else {
			map.at(key) = true;
			return true;
		}
	}
}

bool iterate_remap_check(Quad *qs, uint64_t n, bool is_rect, uint32_t v_offset,
                         uint32_t wing_width, uint32_t q_side_len, edge_list_map &map) {
// iterate over all the edges in the quadrant, remapping the vertex ids to their original values

for (uint64_t i = 0; i < n; ++i) {
		Quad &q = qs[i];
		for (uint32_t j = 0; j < q.nnz; ++j) {
			uint32_t src = q.edges[j * 2];
			uint32_t dest = q.edges[j * 2 + 1];
			if (is_rect) {
//				 tail
				uint32_t u = q.qx + src;
				uint32_t v = q.qy + dest + v_offset;
				bool res = check_edge(u, v, map);
				if (not res) {
					fmt::print("src, dest {}: {}\n", src, dest );
					fmt::print("q.qx, q.qy, q_idx: {} {} {}\n", q.qx, q.qy, q.q_idx);
					fmt::print("u, v: {} {}\n", u, v);

					return false;
				}
			} else {
				// wings
				uint32_t u = q_side_len * q.qx + src;
				uint32_t v = q_side_len * q.qy + dest + v_offset;
				bool res = check_edge(u, v, map);
				if (not res) { return false; }
			}
		}
	}
	return true;
}

bool check_all_edges(Quad *left_wing_qs, uint32_t n_nnz_lw_qs,
                     Quad *right_wing_qs, uint32_t n_nnz_rw_qs,
                     Quad *tail_rects, uint32_t n_nnz_tail_qs,
                     std::vector<std::pair<uint32_t, uint32_t>> edge_list,
                     uint32_t wing_width, uint32_t q_side_len) {

	// construct a map from each to a boolean value
	edge_list_map seen;
	for (uint64_t i = 0; i < edge_list.size(); ++i) {
		seen.insert({{edge_list[i].first, edge_list[i].second},
		             false});
	}
	assert (seen.size() == edge_list.size());
	bool res = false;
	// check that each edge is seen exactly once
	fmt::print("Checking left wing quads..\n");
	if (!iterate_remap_check(left_wing_qs, n_nnz_lw_qs, false, 0, wing_width, q_side_len, seen)) return false;
	fmt::print("Checking right wing quads..\n");
	if (!iterate_remap_check(right_wing_qs, n_nnz_rw_qs, false, wing_width, wing_width, q_side_len, seen)) return false;
	fmt::print("Checking tail rectangles..\n");
	if (!iterate_remap_check(tail_rects, n_nnz_tail_qs, true, 0, wing_width, q_side_len, seen)) return false;

	fmt::print("Checking all edges seen..\n");
	// finally, iterate over the edge map and check that all edges have been seen
	for (const auto &kv: seen) {
		if (!kv.second) {
			fmt::print("kv: {}\n", kv);
			return false;
		}
	}
	return true;
}

/**
 * Preprocesses a graph that has been reordered using the Parallel-Slashburn vertex ordering
 * Preprocessing includes:
 * Splitting the graph into Left-Wing, Right-Wing, and Tail
 * Each of the three sections will be split into quadrants
 * Quadrants within a section will be traversed using a hilbert curve
 * Edges within a quadrant will be traversed using a hilbert curve
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
	opterr = 0;
	int opt;
	uint32_t n = 0; // n vertices
	uint64_t m = 0; // n edges
	int num_expts = 0;
//	int n_threads = 0;
	std::string graph_name = "";
	std::string data_dir = "";
	uint32_t q_side_len = 0;
	uint32_t wing_width = 0;
	bool debug = false;

	while ((opt = getopt(argc, argv, "bg:d:l:e:w:")) != -1) {
		switch (opt) {
			case 'b':
				debug = !debug;
				break;
			case 'g':
				graph_name = optarg;
				break;
			case 'd':
				data_dir = optarg;
				break;
			case 'l':
				q_side_len = atoi(optarg);
			case 'w':
				wing_width = atoi(optarg);
//			case 't':
//				n_threads = atoi(optarg);
//			case 'e':
//				num_expts = atoi(optarg);
		}
	}
	assert(is_power_of_2(q_side_len));
	// read the original edgelist of the graph
	std::vector<uint32_t> iso_map;
	std::vector<std::pair<uint32_t, uint32_t>> edge_list;

	std::string graphs_dir = fmt::format("{}/{}", data_dir, "graphs");
	std::string graph_dir = fmt::format("{}/{}", graphs_dir, graph_name);
	std::string binary_order_path = fmt::format("{}/{}", graph_dir, "parsb.bin");
	std::string binary_edge_list_path = fmt::format("{}/{}", graph_dir, "comp.bin");
	read_binary_container<std::vector<uint32_t>>(binary_order_path, iso_map);
	read_binary_container<std::vector<std::pair<uint32_t, uint32_t>>>(binary_edge_list_path, edge_list);

	n = iso_map.size();
	m = edge_list.size();

	// for debug, maintain an edge set
	std::set<std::pair<uint32_t, uint32_t>> edge_set;
	fmt::print("n: {}\n", n);
	fmt::print("m: {}\n", m);
	// remap the edges of the graph and sort by src, dest
//#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < m; ++i) {
		uint32_t src = edge_list[i].first;
		uint32_t dest = edge_list[i].second;

		edge_list[i].first = iso_map[src];
		edge_list[i].second = iso_map[dest];
		edge_set.insert({iso_map[src], iso_map[dest]});
	}



	// the slashburn isomorphism map, translate the edgelist, and sort by src, dest
//	std::sort(dpl::execution::par_unseq, edge_list.begin(), edge_list.end());
	ips4o::parallel::sort(edge_list.begin(), edge_list.end());

	// read the indegrees, outdegrees of the vertices of the graph, and translate using the
	// slashburn isomap
	std::string out_deg_path = fmt::format("{}/{}", graph_dir, "out_degs");
	std::string in_deg_path = fmt::format("{}/{}", graph_dir, "in_degs");

	std::vector<uint32_t> out_degs(n);
	std::vector<uint32_t> in_degs(n);

	read_degs(out_deg_path, in_deg_path, out_degs, in_degs, n, iso_map);
//	fmt::print("out_degs: {}\n", out_degs);
//	fmt::print("in_degs: {}\n", in_degs);

	// create the in-csr and out-csr representation
	CSR in_csr = CSR(n, m, true);
	CSR out_csr = CSR(n, m, false);

	out_csr.par_populate(out_degs, edge_list);



	// before populating in the in-csr, sort the edgelist by ascending (dest, src)
	ips4o::parallel::sort(
		edge_list.begin(),
		edge_list.end(),
		[](std::pair<uint32_t, uint32_t> l, std::pair<uint32_t, uint32_t> r) -> bool {
			if (l.second == r.second) { return l.first < r.first; }
			return l.second < r.second;
		}
	);
	in_csr.par_populate(in_degs, edge_list);
	uint64_t n_edges_lwing = 0;
	uint64_t n_edges_rwing = 0;
	uint64_t n_edges_tail = 0;

	// separate the edges into left wing, right wing, and tail; first count the number of
	// edges in each section
	// left wing
#pragma omp parallel for schedule(static) reduction(+:n_edges_lwing)
	for (uint32_t u = 0; u < n; ++u) {
		uint64_t start = out_csr.index[u];
		uint64_t end = out_csr.index[u + 1];
		for (uint64_t offset = start; offset < end; ++offset) {
			uint32_t v = out_csr.neighbours[offset];
			if (v >= wing_width) { break; }
			++n_edges_lwing;
		}
	}

	// right wing
#pragma omp parallel for schedule(static) reduction(+:n_edges_rwing)
	for (uint32_t u = 0; u < wing_width; ++u) {
		uint64_t start = out_csr.index[u];
		uint64_t end = out_csr.index[u + 1];
		// only consider the in_edges whose source id is greater than the wing_width
		uint64_t offset = lower_bound<uint32_t, uint64_t>(
			wing_width,
			out_csr.neighbours,
			start,
			end
		);
		for (uint64_t i = offset; i < end; ++i) {
			uint32_t v = out_csr.neighbours[i];
			++n_edges_rwing;
		}
	}
	// tail
#pragma omp parallel for schedule(static) reduction(+:n_edges_tail)
	for (uint32_t u = wing_width; u < n; ++u) {
		uint64_t start = out_csr.index[u];
		uint64_t end = out_csr.index[u + 1];
		// only consider the in_edges whose source id is greater than the wing_width
		uint64_t offset = lower_bound<uint32_t, uint64_t>(
			wing_width,
			out_csr.neighbours,
			start,
			end
		);
		for (uint64_t i = offset; i < end; ++i) {
			uint32_t v = out_csr.neighbours[i];
			++n_edges_tail;
		}
	}
	fmt::print("n_edges_lwing: {}\n", n_edges_lwing);
	fmt::print("n_edges_rwing: {}\n", n_edges_rwing);
	fmt::print("n_edges_tail: {}\n", n_edges_tail);
	fmt::print("q_side_len: {}\n", q_side_len);
	fmt::print("wing_width: {}\n", wing_width);

	// stripe_len divides the wings into separate units of computation
	// (horizontal stripes for the left wing, vertical stripes for the right wing)
	// stripe_len is equal to the largest power of 2 smaller than the wing_width
	uint32_t stripe_len = hyperceiling(wing_width) / 2;
	fmt::print("stripe_len: {}\n", stripe_len);
	// compute the number of quadrants for each section (left, right wing, tail) of the adjacency matrix
	// and compute the number of edges per quadrant (to preallocate the edge array size per quadrant)

	// compute the number of stripes in each wing
	uint32_t n_stripes_in_left_wing = (n + stripe_len - 1) / stripe_len;
	// since we'll compute the edges in the densest quadrant when iterating over the
	// edges of the left wing, we can ignore these when iterating over the right wing
	uint32_t right_wing_width = n - wing_width;
	uint32_t n_stripes_in_right_wing = (right_wing_width + stripe_len - 1) / stripe_len;

	// compute the number of quadrants in the left wing
	uint32_t n_qs_left_wing = 0;
	uint32_t n_qs_per_side = (n + q_side_len - 1) / q_side_len;
	std::vector<uint32_t> n_qs_in_left_wing_row(n_qs_per_side);
	fmt::print("n_qs_per_side: {}\n", n_qs_per_side);

#pragma omp parallel for schedule(static, 1) reduction(+:n_qs_left_wing)
	for (uint32_t row_idx = 0; row_idx < n; row_idx += q_side_len) {
		uint32_t max_dest_in_row = 0;
		for (uint32_t u = row_idx; u < row_idx + q_side_len && u < n; ++u) {
			uint64_t start = out_csr.index[u];
			uint64_t end = out_csr.index[u + 1];

			// if the vertex does not have any neighbours that reside in the left wing,
			// continue
			if (out_csr.neighbours[start] >= wing_width || start == end) { continue; }

			// find the index of the largest destination whose id is lesser than the wing width
			uint64_t offset = lower_bound<uint32_t, uint64_t>(
				wing_width,
				out_csr.neighbours,
				start,
				end
			);
			uint32_t largest_dest_seen = out_csr.neighbours[offset - 1];
			if (largest_dest_seen > max_dest_in_row) max_dest_in_row = largest_dest_seen;
		}


// TODO is this offset below mandatoyr
		uint32_t n_quads_in_row = (max_dest_in_row + q_side_len) / q_side_len;
		n_qs_in_left_wing_row[row_idx / q_side_len] = n_quads_in_row;
		n_qs_left_wing += n_quads_in_row;
	}
	fmt::print("n_qs_left_wing: {}\n", n_qs_left_wing);
	std::vector<uint32_t> cumulative_n_qs_per_row_in_left_wing(n_qs_per_side + 1);
	fmt::print("n_qs_in_left_wing_row: {}\n", n_qs_in_left_wing_row);
	// now that we know the number of quadrants in the right wing create an array of quadrants
	std::inclusive_scan(
		dpl::execution::par_unseq,
		n_qs_in_left_wing_row.begin(),
		n_qs_in_left_wing_row.end(),
		cumulative_n_qs_per_row_in_left_wing.begin() + 1,
		[](uint32_t l, uint32_t r) -> uint32_t { return l + r; },
		0
	);
	fmt::print("cumulative_n_qs_per_row_in_left_wing: {}\n", cumulative_n_qs_per_row_in_left_wing);
	Quad *left_wing_qs = new Quad[n_qs_left_wing]();

	// assign the quad coordinates for the quadrants in the left wing
#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < n_qs_per_side; i++) {
		uint32_t n_qs = n_qs_in_left_wing_row[i];
		for (uint32_t j = 0; j < n_qs; ++j) {
			uint32_t quad_idx = cumulative_n_qs_per_row_in_left_wing[i] + j;
			Quad &q = left_wing_qs[quad_idx];
			q.q_idx = quad_idx;
			q.qx = i;
			q.qy = j;
		}
	}

//	std::vector<Quad> left_wing_qs(n_qs_left_wing);
	// compute the number of edges that will reside in each quadrant of the left wing
	// this can be done in parallel since each thread exclusively inspects a quadrant
	// row
#pragma omp parallel for schedule(static) reduction(+:n_qs_left_wing)
	for (uint32_t row_idx = 0; row_idx < n; row_idx += q_side_len) {
		for (uint32_t u = row_idx; u < row_idx + q_side_len && u < n; ++u) {
			uint64_t start = out_csr.index[u];
			uint64_t end = out_csr.index[u + 1];

			// skip vertices that do not have any neighbours that reside in the left wing
			if (out_csr.neighbours[start] >= wing_width || start == end) { continue; }

			uint32_t quad_row = u / q_side_len;
			for (uint64_t offset = start; offset < end; offset++) {
				uint32_t v = out_csr.neighbours[offset];
				if (v >= wing_width) break;
				uint32_t quad_col = v / q_side_len;
				uint32_t quad_idx = cumulative_n_qs_per_row_in_left_wing[quad_row] + quad_col;
				Quad &q = left_wing_qs[quad_idx];
//				q.q_idx = quad_idx;
//				q.qx = quad_row;
//				q.qy = quad_col;
				q.nnz++;
			}
		}
	}

	// now that we know the number of edges in each quadrant, preallocate
	// the quad local edges array to that size times 2 (since edges will be stored
	// flattened)
#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < n_qs_left_wing; ++i) {
		Quad &q = left_wing_qs[i];
		if (q.nnz == 0) {
//			 if this quadrant is empty, it should not exist
//			assert(false);
		}
		q.edges = new uint32_t[q.nnz * 2]();
	}

	// (logically) separate the left wing into rows of size stripe_len
	// within a stripe, each thread will populate a quadrant row by iterating over the
	// edges in that quadrant row
	uint32_t *n_edges_seen_in_left_wing_q = new uint32_t[n_qs_left_wing]();

//	fmt::print("stripe_len: {}\n", stripe_len);
//	fmt::print("n_stripes_per_side: {}\n", n_stripes_in_left_wing);
//	fmt::print("n_stripes_per_side: {}\n", n_stripes_in_right_wing);
	for (uint32_t stripe_idx = 0; stripe_idx < n_stripes_in_left_wing; stripe_idx++) {
		uint32_t stripe_start = stripe_idx * stripe_len;
		uint32_t stripe_end = (stripe_idx + 1) * stripe_len;
//		fmt::print("stripe_Start: {}\n", stripe_start);
//		fmt::print("stripe_end: {}\n", stripe_end);

#pragma omp parallel for schedule(static, 1)
		for (uint32_t row_idx = stripe_start; row_idx < stripe_end; row_idx += q_side_len) {
			for (uint32_t u = row_idx; u < row_idx + q_side_len && u < n; ++u) {
				uint64_t start = out_csr.index[u];
				uint64_t end = out_csr.index[u + 1];

				// skip vertices that do not have any neighbours that reside in the left wing
				if (out_csr.neighbours[start] >= wing_width || start == end) { continue; }

				uint32_t quad_row = u / q_side_len;
				for (uint64_t offset = start; offset < end; offset++) {
					uint32_t v = out_csr.neighbours[offset];
					if (!edge_set.count({u, v})) {
						fmt::print("left: {} {}\n", u, v);
					}
					if (v >= wing_width) break;
					uint32_t quad_col = v / q_side_len;
					uint32_t quad_idx = cumulative_n_qs_per_row_in_left_wing[quad_row] + quad_col;
					Quad &q = left_wing_qs[quad_idx];
					// compute the local (x, y) coordinate given (u, v)
					auto p = edge_offset(u, v, q_side_len);
					uint32_t qu = p.first;
					uint32_t qv = p.second;
					// insert the offset (u, v) to the flattened edges array of q
					uint32_t edge_offset_into_q = n_edges_seen_in_left_wing_q[quad_idx];
					q.edges[edge_offset_into_q] = qu;
					q.edges[edge_offset_into_q + 1] = qv;
					n_edges_seen_in_left_wing_q[quad_idx] += 2;
				}
			}
		}
	}



	// compute the number of quadrants in the right wing
	uint32_t n_qs_right_wing = 0;
	uint32_t n_qs_in_right_wing_side = (wing_width + q_side_len - 1) / q_side_len;
	std::vector<uint32_t> n_qs_in_right_wing_row(n_qs_in_right_wing_side);
	fmt::print("n_qs_in_right_wing_side: {}\n", n_qs_in_right_wing_side);
	fmt::print("n_qs_in_right_wing_row: {}\n", n_qs_in_right_wing_row);

#pragma omp parallel for schedule(static, 1) reduction(+:n_qs_right_wing)
	for (uint32_t row_idx = 0; row_idx < wing_width; row_idx += q_side_len) {
		uint32_t max_dest_in_row = 0;
		for (uint32_t u = row_idx; u < row_idx + q_side_len && u < wing_width; ++u) {
			uint64_t start = out_csr.index[u];
			uint64_t end = out_csr.index[u + 1];

			if (start == end) continue;

			// get the largest id neighbour of u
			uint32_t largest_dest_seen = out_csr.neighbours[end - 1];
			if (u == 1026) {
				fmt::print("largest_dest_seen: {}\n", largest_dest_seen);
			}
			// if the vertex does not have any neighbours that reside in the right wing,
			// continue
			if (largest_dest_seen < wing_width) continue;
			largest_dest_seen -= wing_width; // offset the neighbour id
			if (largest_dest_seen > max_dest_in_row) max_dest_in_row = largest_dest_seen;
		}
		uint32_t n_quads_in_row = (max_dest_in_row + q_side_len) / q_side_len;
		n_qs_in_right_wing_row[row_idx / q_side_len] = n_quads_in_row;
		n_qs_right_wing += n_quads_in_row;
	}
	fmt::print("n_qs_right_wing: {}\n", n_qs_right_wing);
	std::vector<uint32_t> cumulative_n_qs_per_row_in_right_wing(n_qs_in_right_wing_side + 1);

	fmt::print("n_qs_in_right_wing_row: {}\n", n_qs_in_right_wing_row);
	// now that we know the number of quadrants in the right wing, create an array of quadrants
	std::inclusive_scan(
		dpl::execution::par_unseq,
		n_qs_in_right_wing_row.begin(),
		n_qs_in_right_wing_row.end(),
		cumulative_n_qs_per_row_in_right_wing.begin() + 1,
		[](uint32_t l, uint32_t r) -> uint32_t { return l + r; },
		0
	);
	fmt::print("cumulative_n_qs_per_row_in_right_wing: {}\n", cumulative_n_qs_per_row_in_right_wing);
	Quad *right_wing_qs = new Quad[n_qs_right_wing]();

	// assign the quad coordinates for the quadrants in the right wing
#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < n_qs_in_right_wing_side; i++) {
		uint32_t n_qs = n_qs_in_right_wing_row[i];
		for (uint32_t j = 0; j < n_qs; ++j) {
			uint32_t quad_idx = cumulative_n_qs_per_row_in_right_wing[i] + j;
//			fmt::print("i, j: {} {} {}\n", i, j, quad_idx);

			Quad &q = right_wing_qs[quad_idx];
			q.q_idx = quad_idx;
			q.qx = i;
			q.qy = j;
		}
	}

	// compute the number of edges that will reside in each quadrant of the right wing
	// this can be done in parallel since each thread will exclusively modify a quadrant
	// row
#pragma omp parallel for schedule(static, 1)
	for (uint32_t row_idx = 0; row_idx < wing_width; row_idx += q_side_len) {
		for (uint32_t u = row_idx; u < row_idx + q_side_len && u < wing_width; ++u) {
			uint64_t start = out_csr.index[u];
			uint64_t end = out_csr.index[u + 1];

			// skip vertices that do not have any neighbours that reside in the right wing
			if (out_csr.neighbours[end - 1] < wing_width || start == end) { continue; }

			// find the index of the largest destination whose id is lesser than the wing width
			uint64_t edge_start = lower_bound<uint32_t, uint64_t>(
				wing_width,
				out_csr.neighbours,
				start,
				end
			);


			uint32_t quad_row = u / q_side_len;
			for (uint64_t offset = edge_start; offset < end; offset++) {
				uint32_t v = out_csr.neighbours[offset];
				uint32_t quad_col = (v - wing_width) / q_side_len;
				uint32_t quad_idx = cumulative_n_qs_per_row_in_right_wing[quad_row] + quad_col;
				Quad &q = right_wing_qs[quad_idx];
//				q.q_idx = quad_idx;
//				q.qx = quad_row;
//				q.qy = quad_col;
				q.nnz++;
			}
		}
	}

	// now that we know the number of edges in each quadrant, preallocate
	// the quad local edges array to that size times 2 (since edges will be stored
	// flattened)
#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < n_qs_right_wing; ++i) {
		Quad &q = right_wing_qs[i];
		if (q.nnz == 0) {
//			 if this quadrant is empty, it should not exist
//			assert(false);
		}
//		fmt::print("{} {} {} {}\n", q.q_idx,q.qx, q.qy, q.nnz);
		q.edges = new uint32_t[q.nnz * 2]();
	}

	//	right 1282 2940
	uint32_t t1 = 1282;
	uint64_t o1 = out_csr.index[t1];
	uint64_t o2 = out_csr.index[t1 + 1];
	for (uint64_t o3 = o1; o3 < o2; o3++) {
		fmt::print("out_csr.neighbours[o3]: {}\n", out_csr.neighbours[o3]);
	}



//	 (logically) separate the right wing into columns of size stripe_len
//	 within a stripe, each thread will populate a quadrant row by iterating over the
//	 edges in that quadrant column
	uint32_t *n_edges_seen_in_right_wing_q = new uint32_t[n_qs_right_wing]();

	for (uint32_t stripe_idx = 0; stripe_idx < n_stripes_in_right_wing; stripe_idx++) {
		uint32_t stripe_start = stripe_idx * stripe_len + wing_width;
		uint32_t stripe_end = (stripe_idx + 1) * stripe_len + wing_width;
//		fmt::print("stripe_Start: {}\n", stripe_start);
//		fmt::print("stripe_end: {}\n", stripe_end);

#pragma omp parallel for schedule(static, 1)
		for (uint32_t col_idx = stripe_start; col_idx < stripe_end; col_idx += q_side_len) {
			for (uint32_t v = col_idx; v < col_idx + q_side_len && v < n; ++v) {
				uint64_t start = in_csr.index[v];
				uint64_t end = in_csr.index[v + 1];

				// skip vertices that do not have any in-neighbours that reside in the right wing
				if (in_csr.neighbours[start] >= wing_width || start == end) { continue; }

				uint32_t quad_col = (v - wing_width) / q_side_len;
				for (uint64_t offset = start; offset < end; offset++) {
					uint32_t u = in_csr.neighbours[offset];
					if (u >= wing_width) break;
					uint32_t quad_row = u / q_side_len;
					uint32_t quad_idx = cumulative_n_qs_per_row_in_right_wing[quad_row] + quad_col;
					Quad &q = right_wing_qs[quad_idx];


					// compute the local (x, y) coordinate given (u, v)

					auto p = edge_offset(u, v - wing_width, q_side_len);
					uint32_t qu = p.first;
					uint32_t qv = p.second;
					// insert the offset (u, v) to the flattened edges array of q
					uint32_t edge_offset_into_q = n_edges_seen_in_right_wing_q[quad_idx];
//					assert (edge_offset_into_q < q.nnz * 2);
					q.edges[edge_offset_into_q] = qu;
					q.edges[edge_offset_into_q + 1] = qv;
					n_edges_seen_in_right_wing_q[quad_idx] += 2;
				}
			}
		}
	}

/**
	 * this solution takes the cumulative sum of number of edges
	 * for each row in the tail, and divides up the total number of edges among the given
	 * number of quadrants
	 */
	uint32_t tail_size = n - wing_width;
	uint32_t *tail_degrees = new uint32_t[tail_size]();
	uint64_t total_edges_in_tail = 0;
	// compute the number of edges for each row in the tail
#pragma omp parallel for schedule(static) reduction(+:total_edges_in_tail)
	for (uint32_t u = wing_width; u < n; ++u) {
		uint64_t start = out_csr.index[u];
		uint64_t end = out_csr.index[u + 1];
		if (start == end || end < start) continue;
		uint64_t offset = lower_bound<uint32_t, uint64_t>(
			wing_width,
			out_csr.neighbours,
			start,
			end
		);
		if (offset == end) continue;
		for (uint32_t v = offset; v < end; ++v) {
			++total_edges_in_tail;
			++tail_degrees[u - wing_width];
		}
	}
//	(max_dest_in_row + q_side_len - 1) / q_side_len;
	uint32_t n_rects_in_tail = (tail_size + q_side_len - 1) / q_side_len;
	uint32_t *rect_lbs = new uint32_t[n_rects_in_tail]();
	uint32_t *rect_ubs = new uint32_t[n_rects_in_tail]();
	Quad *tail_rects = new Quad[n_rects_in_tail]();

	// split the tail into rectangles of width q_side_len
	// if edges exist outside the bounds of the rectangle, increase its length appropriately
	for (uint32_t q = wing_width; q < n; q += q_side_len) {

		uint32_t min_in_rect = -1;
		uint32_t max_in_rect = 0;
		uint32_t q_idx = (q - wing_width) / q_side_len;
		uint32_t q_ub = q + q_side_len >= n ? n : q + q_side_len;
#pragma omp parallel for schedule(static) reduction (max:max_in_rect) reduction(min:min_in_rect)
		for (uint32_t u = q; u < q_ub; ++u) {
			// get the x-bounds of this rectangle
			uint32_t x_start = in_csr.index[u];
			uint32_t x_end = in_csr.index[u + 1];

			uint32_t quad_lb = q;
			uint32_t quad_ub = q + q_side_len;
			if (quad_ub >= n) quad_ub = n;
			uint32_t max_in_neigh = in_csr.neighbours[x_end - 1];
			// find any tail in-neighbours that may be out of bounds from above
			uint32_t min_in_neigh = get_min_id_in_range(u, wing_width, quad_lb, in_csr);
			if (max_in_neigh > quad_ub)
				if (max_in_neigh > max_in_rect) max_in_rect = max_in_neigh;
			if (min_in_neigh < quad_lb && min_in_neigh >= wing_width)
				if (min_in_neigh < min_in_rect) min_in_rect = min_in_neigh;
		}
		rect_lbs[q_idx] = min_in_rect;
		rect_ubs[q_idx] = max_in_rect;
	}

	// if no values have been assigned to the bounds assign them to the original quad bounds
#pragma omp parallel for schedule(static)
	for (uint32_t q = wing_width; q < n; q += q_side_len) {
		uint32_t q_idx = (q - wing_width) / q_side_len;
		if (rect_lbs[q_idx] == -1) rect_lbs[q_idx] = q;
		uint32_t ub = q + q_side_len > n ? n : q + q_side_len;
		if (rect_ubs[q_idx] == 0) rect_ubs[q_idx] = ub;
	}
	rect_lbs[0] = wing_width;
	rect_ubs[n_rects_in_tail - 1] = n - 1;

	// assign the quad coordinates for the quadrants in the tail

	// compute the number of edges in each tail rectangles to preallocate their edges arrays
	// count the number of edges that lie within each tail rectangle
	// important - tail upper and lower bounds for rectangle are inclusive (length bounds)
	// the q_idx of all tail rectangles can be omitted, since the order of iterations of
	// rects is defined by their qy coordinate
	// as such, store the bounds of the rectangle's length in qx, q_idx
	// i.e. a tail rectangle, r, is defined by its anchor point (r.qx, r.qy)
	// r's width == q_side_len
	// r's length == r.q_idx - r.qx
	uint32_t total = 0;
	uint32_t t_total = 0;
#pragma omp parallel for schedule(static) reduction(+:total)
	for (uint32_t q = wing_width; q < n; q += q_side_len) {
//		fmt::print("q: {}\n", q);
		uint32_t q_idx = (q - wing_width) / q_side_len;
		uint32_t q_ub = q + q_side_len >= n ? n : q + q_side_len;
		uint32_t n_neighbours_in_rect = 0;
		uint32_t x_lb = rect_lbs[q_idx];
		uint32_t x_ub = rect_ubs[q_idx];

		for (uint32_t u = x_lb; u <= x_ub; ++u) {

			uint32_t n_neighbours_in_row = n_out_neighbours_between(u, q, q_ub, out_csr);
			n_neighbours_in_rect += n_neighbours_in_row;
			uint64_t start = out_csr.index[u];
			uint64_t end = out_csr.index[u + 1];
		}
		Quad &r = tail_rects[q_idx];
		r.nnz = n_neighbours_in_rect;
		r.qx = x_lb;
		r.q_idx = x_ub + 1;
		r.qy = q;
		total += n_neighbours_in_rect;
	}

	// preallocate each edges array in the tail's rectangles
#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < n_rects_in_tail; ++i) {
		Quad &q = tail_rects[i];
		if (q.nnz == 0) {
//			 if this quadrant is empty, it should not exist
//			assert(false);
		}
//		fmt::print("{} {} {} {}\n", q.q_idx,q.qx, q.qy, q.nnz);
		q.edges = new uint32_t[q.nnz * 2]();
	}

	// assign all edges within a quadrant to its corresponding edges array
	uint32_t *n_edges_seen_in_tail_rect = new uint32_t[n_rects_in_tail]();
	uint64_t test = 0;
#pragma omp parallel for schedule(static) reduction(+:test)
	for (uint32_t quad_idx = 0; quad_idx < n_rects_in_tail; ++quad_idx) {
		Quad &q = tail_rects[quad_idx];
		// care has to be taken to compute the offsets of the rectangle - it is _not_ computed
		// q_side_len, as with quadrants in the left, right wings
		// instead the offset is computed using the largest power of 2 that is greater than
		// the rectangle's length

		uint32_t rect_length = q.q_idx - q.qx;
		uint32_t rect_offset = hyperceiling(rect_length);
		for (uint32_t u = q.qx; u <= q.q_idx; ++u) {
			uint32_t start = out_csr.index[u];
			uint32_t end = out_csr.index[u + 1];

			// find the first out neighbour of u in this rectangle
			uint64_t offset = lower_bound<uint32_t, uint64_t>(
				q.qy,
				out_csr.neighbours,
				start,
				end
			);
			if (offset == end) continue;

			for (uint64_t out_edge_offset = offset; out_edge_offset < end; ++out_edge_offset) {
				uint32_t v = out_csr.neighbours[out_edge_offset];
				if (!edge_set.count({u, v})) {
					fmt::print("tail: {} {}\n", u, v);
				}
				if (v >= q.qy + q_side_len) break;
				auto p = edge_offset(u - q.qx, v - q.qy, rect_offset);
				uint32_t qu = p.first;
				uint32_t qv = p.second;
//				insert the offset (u, v) to the flattened edges array of q
				uint32_t edge_offset_into_q = n_edges_seen_in_tail_rect[quad_idx];
				q.edges[edge_offset_into_q] = qu;
				q.edges[edge_offset_into_q + 1] = qv;
				n_edges_seen_in_tail_rect[quad_idx] += 2;
				test += 1;
			}
		}
	}
	// assert that all edges in left, right wings and tail have been assigned to a quadrant
	uint64_t lw_count_test = 0;
	uint64_t rw_count_test = 0;
	uint64_t tail_count_test = 0;

#pragma omp parallel for schedule(static) reduction(+:lw_count_test)
	for (uint32_t i = 0; i < n_qs_left_wing; ++i)
		lw_count_test += left_wing_qs[i].nnz;

#pragma omp parallel for schedule(static) reduction(+:rw_count_test)
	for (uint32_t i = 0; i < n_qs_right_wing; ++i)
		rw_count_test += right_wing_qs[i].nnz;

#pragma omp parallel for schedule(static) reduction(+:tail_count_test)
	for (uint32_t i = 0; i < n_rects_in_tail; ++i)
		tail_count_test += tail_rects[i].nnz;

	assert(lw_count_test == n_edges_lwing);
	assert(rw_count_test == n_edges_rwing);
	assert(tail_count_test == n_edges_tail);
	assert(lw_count_test + rw_count_test + tail_count_test == m);

	// hilbert order the edges within all quadrants
#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < n_qs_left_wing; ++i)
		horder_edges_in_q(left_wing_qs[i], q_side_len);

#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < n_qs_right_wing; ++i)
		horder_edges_in_q(right_wing_qs[i], q_side_len);

#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < n_rects_in_tail; ++i) {
		Quad &q = tail_rects[i];
		uint32_t rect_length = q.q_idx - q.qx;
		uint32_t rect_offset = hyperceiling(rect_length);
		horder_edges_in_q(q, rect_offset);
	}


// hilbert order the quadrants in the wings
// quadrants are chunked into rows of size stripe_len (the largest power of 2 that is
// smaller than the wing width)

	// count the number of quadrants within each stripe
	uint32_t n_quads_per_stripe = (stripe_len + q_side_len - 1) / q_side_len;
	uint32_t *cumulative_n_quads_per_left_wing_stripe = new uint32_t[n_stripes_in_left_wing + 1]();

	uint32_t left_wing_stripe_idx = 1;
	for (uint32_t i = n_quads_per_stripe; i < n_qs_per_side; i += n_quads_per_stripe) {
//		fmt::print("i: {}\n", i);
		cumulative_n_quads_per_left_wing_stripe[left_wing_stripe_idx] = cumulative_n_qs_per_row_in_left_wing[i];
		++left_wing_stripe_idx;
	}
	cumulative_n_quads_per_left_wing_stripe[n_stripes_in_left_wing] = cumulative_n_qs_per_row_in_left_wing[n_qs_per_side];
//	print_arr<uint32_t>(cumulative_n_quads_per_left_wing_stripe, n_stripes_in_left_wing + 1);

	fmt::print("n_stripes_in_left_wing: {}\n", n_stripes_in_left_wing);
	fmt::print("n_quads_per_stripe: {}\n", n_quads_per_stripe);
	// hilbert sort all quadrants within a stripe
	for (uint32_t i = 0; i < n_stripes_in_left_wing; ++i) {
		uint32_t stripe_start = cumulative_n_quads_per_left_wing_stripe[i];
		uint32_t stripe_end = cumulative_n_quads_per_left_wing_stripe[i + 1];
//		fmt::print("stripe_start, stripe_end: {} {}\n", stripe_start, stripe_end);
		// assign the hilbert idcs
#pragma omp parallel for schedule(static)
		for (uint32_t j = stripe_start; j < stripe_end; ++j) {
			Quad &q = left_wing_qs[j];

			// offset the stripe
			uint32_t qx = q.qx - (n_quads_per_stripe * i);

			// the coordinate system used to compute the hilbert index of each offset stripe
			// is defined using the number of quads per stripe side length *2
			uint32_t hilbert_side_len = n_quads_per_stripe * 2;
			uint32_t h_idx = xy2d(hilbert_side_len, qx, q.qy);
			q.q_idx = h_idx;
		}

		// sort the quadrants within this horizontal stripe by ascending hilbert index
		std::sort(
			dpl::execution::par_unseq,
			left_wing_qs + stripe_start,
			left_wing_qs + stripe_end,
			[](const Quad &l, const Quad &r) { return l.q_idx < r.q_idx; }
		);
	}

	// before reordering the quadrants in the right wing, sort them by ascending qy
	// and cumulatively sum:
	// - the number of quadrants in each quadrant column
	// - the number of quadrants in each vertical stripes
	std::sort(
		dpl::execution::par_unseq,
		right_wing_qs,
		right_wing_qs + n_qs_right_wing,
		[](const Quad &l, const Quad &r) {
			if (l.qy == r.qy) return l.qx < r.qx;
			return l.qy < r.qy;
		}
	);
	// get the number of quadrants in each quadrant column
	uint32_t n_cols_in_right_wing = (right_wing_width + q_side_len - 1) / q_side_len;
//	uint32_t *n_quads_in_right_wing_col = new uint32_t[n_cols_in_right_wing]();
	std::vector<uint32_t> n_quads_in_right_wing_col(n_cols_in_right_wing);
	uint32_t curr_col = 0;
	uint32_t quad_count = 0;
	for (uint32_t j = 0; j < n_qs_right_wing; ++j) {
		Quad &q = right_wing_qs[j];
		if (q.qy != curr_col) {
			n_quads_in_right_wing_col[curr_col] = quad_count;
			curr_col = q.qy;
			quad_count = 0;
		}
		quad_count++;
	}
	n_quads_in_right_wing_col[n_cols_in_right_wing - 1] = quad_count;
	std::vector<uint32_t> cumulative_n_qs_per_col_in_right_wing(n_cols_in_right_wing + 1);
//	print_arr<uint32_t>(n_quads_in_right_wing_col, n_cols_in_right_wing);
	fmt::print("n_quads_in_right_wing_col: {}\n", n_quads_in_right_wing_col);
	fmt::print("n_cols_in_right_wing: {}\n", n_cols_in_right_wing);
	std::inclusive_scan(
		dpl::execution::par_unseq,
		n_quads_in_right_wing_col.begin(),
		n_quads_in_right_wing_col.end(),
		cumulative_n_qs_per_col_in_right_wing.begin() + 1,
		[](uint32_t l, uint32_t r) -> uint32_t { return l + r; },
		0
	);
	fmt::print("cumulative_n_qs_per_col_in_right_wing: {}\n", cumulative_n_qs_per_col_in_right_wing);
	fmt::print("n_quads_per_stripe: {}\n", n_quads_per_stripe);

	uint32_t *cumulative_n_quads_per_right_wing_stripe = new uint32_t[n_stripes_in_right_wing + 1]();

	uint32_t n_qs_in_right_wing_length = ((n - wing_width) + q_side_len - 1) / q_side_len;
	uint32_t right_wing_stripe_idx = 1;

	for (uint32_t i = n_quads_per_stripe; i < n_qs_in_right_wing_length; i += n_quads_per_stripe) {
//		fmt::print("i: {}\n", i);
		cumulative_n_quads_per_right_wing_stripe[right_wing_stripe_idx] = cumulative_n_qs_per_col_in_right_wing[i];
		++right_wing_stripe_idx;
	}

	cumulative_n_quads_per_right_wing_stripe[n_stripes_in_right_wing] = cumulative_n_qs_per_col_in_right_wing[n_qs_in_right_wing_length];
//	print_arr<uint32_t>(cumulative_n_quads_per_right_wing_stripe, n_stripes_in_right_wing + 1);

	// hilbert sort all quadrants within a stripe

	for (uint32_t i = 0; i < n_stripes_in_right_wing; ++i) {
		uint32_t stripe_start = cumulative_n_quads_per_right_wing_stripe[i];
		uint32_t stripe_end = cumulative_n_quads_per_right_wing_stripe[i + 1];
//		fmt::print("stripe_start, stripe_end: {} {}\n", stripe_start, stripe_end);
		// assign the hilbert idcs
#pragma omp parallel for schedule(static)
		for (uint32_t j = stripe_start; j < stripe_end; ++j) {
			Quad &q = right_wing_qs[j];
			// offset the vertical stripe
			uint32_t qy = q.qy - (n_quads_per_stripe * i);

			// the coordinate system used to compute the hilbert index of each offset stripe
			// is defined using the number of quads per stripe side length *2
			uint32_t hilbert_side_len = n_quads_per_stripe * 2;

			// since we want right wing stripes to be traversed vertically, swap qx and qy
			// when computing the hilbert index
			uint32_t h_idx = xy2d(hilbert_side_len, qy, q.qx);
			q.q_idx = h_idx;
		}
//		// sort the quadrants within this horizontal stripe by ascending hilbert index
		std::sort(
			dpl::execution::par_unseq,
			right_wing_qs + stripe_start,
			right_wing_qs + stripe_end,
			[](const Quad &l, const Quad &r) { return l.q_idx < r.q_idx; }
		);
	}


	// before writing quadrants to disk, filter out any empty quads
	// filter out fully zero quads
	uint32_t n_nnz_lw_qs = filter_empty_quads(left_wing_qs, n_qs_left_wing);
	uint32_t n_nnz_rw_qs = filter_empty_quads(right_wing_qs, n_qs_right_wing);
	uint32_t n_nnz_tail_qs = filter_empty_quads(tail_rects, n_rects_in_tail);

	// verify all edges assigned
	if (debug) {
		fmt::print("Verifying all edges have been correctly assigned to quadrants..\n");
		bool all_edges_assigned = check_all_edges(
			left_wing_qs, n_nnz_lw_qs,
			right_wing_qs, n_nnz_rw_qs,
			tail_rects, n_nnz_tail_qs,
			edge_list,
			wing_width, q_side_len
		);
		assert(all_edges_assigned);
	} else {
		fmt::print("Debug flag not passed; Skipping verification of edges.\n");
	}


	// write the quadrants to separate binary files
	std::string lw_path = fmt::format("{}/{}", graph_dir, "lw.bin");
	std::string rw_path = fmt::format("{}/{}", graph_dir, "rw.bin");
	std::string tail_path = fmt::format("{}/{}", graph_dir, "tail.bin");
	write_quad_array(lw_path, left_wing_qs, n_nnz_lw_qs, q_side_len, wing_width, n, m, false);
	write_quad_array(rw_path, right_wing_qs, n_nnz_rw_qs, q_side_len, wing_width, n, m, false);
	write_quad_array(tail_path, tail_rects, n_nnz_tail_qs, q_side_len, wing_width, n, m, true);

//	fmt::print("{:*^30}\n", "left wing");
//	for (uint32_t j = 0; j < n_qs_left_wing; ++j) {
//		Quad &q = left_wing_qs[j];
//		fmt::print(" {} {} {} {} \n", q.qx, q.qy, q.q_idx, q.nnz);
//	}
//
//	fmt::print("{:*^30}\n", "right wing");
//
//	for (uint32_t j = 0; j < n_qs_right_wing; ++j) {
//		Quad &q = right_wing_qs[j];
//		fmt::print(" {} {} {} {} \n", q.qx, q.qy, q.q_idx, q.nnz);
//	}
//
//	fmt::print("{:*^30}\n", "tail");
//
//	for (uint32_t j = 0; j < n_rects_in_tail; ++j) {
//		Quad &q = tail_rects[j];
//		fmt::print(" {} {} {} {} \n", q.qx, q.qy, q.q_idx, q.nnz);
//
//	}

	// clean up
//	delete[]left_wing_qs;
	delete[]n_edges_seen_in_left_wing_q;
	delete[]right_wing_qs;
	delete[]n_edges_seen_in_right_wing_q;
	delete[]tail_degrees;
	delete[]rect_lbs;
	delete[]rect_ubs;
	delete[]tail_rects;
	delete[]n_edges_seen_in_tail_rect;
	delete[]cumulative_n_quads_per_left_wing_stripe;
//	delete[]splits;
//	delete[]min_tail_bounds;
//	delete[]max_tail_bounds;
//	delete[]tail_qs;
}
/**
 * // compute the even split points
	uint32_t accum = 0;
	int split_count = 0;
	uint32_t *splits = new uint32_t[n_threads + 1]();
	for (uint32_t u = 0; u < tail_size; ++u) {
		if (accum >= n_edges_per_thread) {
			accum = 0;
			splits[split_count + 1] = u;
			split_count++;
		}
		accum += tail_degrees[u];
	}

	splits[n_threads] = n - wing_width;
//	print_arr<uint32_t>(splits, n_threads + 1);

	// these are the y bounds for each tail rectangle
	// both min and max are needed since rectangles may overlap (i.e. the start point
	// of consecutive rects are _not_ the end point of the previous)
	uint32_t *min_tail_bounds = new uint32_t[n_threads]();
	std::fill_n(min_tail_bounds, n_threads, -1);
	uint32_t *max_tail_bounds = new uint32_t[n_threads]();

	// the splits array defines the boundaries of the square quadrants that will make up the tail
	// however, edges may still exist outside said boundaries
	// compute the length of the rectangles by looking at the min, max ids of the in-neighbours
	// in each rectangle
	for (int i = 0; i < n_threads; ++i) {
		uint32_t quad_start = splits[i];
		uint32_t quad_end = splits[i + 1];
//		fmt::print("quad_start, quad_end {} {}\n", quad_start, quad_end);
		uint32_t mn = -1;
		uint32_t mx = 0;
//#pragma omp parallel for schedule(static) reduction(max:mx) reduction(min:mn)
		for (uint32_t u = quad_start; u < quad_end; ++u) {
			// for each vertex, look at the min, max id of it's in-neighbourhood
			// if this id lies outside of the range of the quadrant, update the bounds of the
			// square to a vertical rectangle
			uint32_t start = in_csr.index[u + wing_width];
			uint32_t end = in_csr.index[u + wing_width + 1];
			uint32_t max_in_neigh = in_csr.neighbours[end - 1];
			uint32_t quad_lower_bound = splits[i] + wing_width;
			uint32_t quad_upper_bound = splits[i + 1] + wing_width;

			// find any tail in-neighbours that may be out of bounds from above
			uint32_t min_in_neigh = get_min_id_in_range(u + wing_width, wing_width + 1, quad_lower_bound, in_csr);

			if (max_in_neigh > quad_upper_bound)
				if (max_in_neigh > mx) mx = max_in_neigh;
			if (min_in_neigh < quad_lower_bound && min_in_neigh > wing_width)
				if (min_in_neigh < mn) mn = min_in_neigh;
		}

		min_tail_bounds[i] = mn;
		max_tail_bounds[i] = mx;
	}
	// offset the start x coord of each rectangle by adding back the wing width
	for (int i = 0; i < n_threads + 1; i++) {
		splits[i] += wing_width;
	}
	// the y lower bound of the first tail rectangle is always the wing width
	min_tail_bounds[0] = wing_width;
	// the y upper bound of the last tail rectangle is always n
	max_tail_bounds[n_threads - 1] = n - 1;

	print_arr<uint32_t>(splits, n_threads + 1);
	print_arr<uint32_t>(min_tail_bounds, n_threads);
	print_arr<uint32_t>(max_tail_bounds, n_threads);

	// count the number of edges that lie within each tail rectangle
	uint32_t total = 0;
//#pragma omp parallel for schedule(static) reduction(+:total)
	for (int t = 0; t < n_threads; t++) {
		uint32_t y_start = splits[t];
		uint32_t y_end = splits[t + 1];
		uint32_t x_start = min_tail_bounds[t];
		uint32_t x_end = max_tail_bounds[t];
		uint32_t n_neighbours_in_rect = 0;
		for (uint32_t u = x_start; u <= x_end; ++u) {
			uint32_t n_neighbours_in_row = n_out_neighbours_between(u, y_start, y_end , out_csr);
			n_neighbours_in_rect += n_neighbours_in_row;
		}
		fmt::print("n_neighbours_in_rect: {}\n", n_neighbours_in_rect);
		total += n_neighbours_in_rect;
	}
	fmt::print("t: {}\n", total);
 */
/**
 * * todo below is an alternative soln to dividing up the work to be done on the tail
	 * this solution divides the tail into quadrants but gets complicated where
	 * neighbouring vertical quadrants may overlap in writes,
 */

// compute the number of quadrants, and number of edges per quadrant in the tail
// compute the number of quadrants in the right wing
//	uint32_t tail_size = n - wing_width;
//	uint32_t n_qs_tail = 0;
//	uint32_t n_qs_in_tail_side = (tail_size + q_side_len - 1) / q_side_len;
//	std::vector<uint32_t> n_qs_in_tail_row(n_qs_in_tail_side);
//	std::vector<uint32_t> tail_lbs(n_qs_in_tail_side); // lower bounds of tail
//	std::vector<uint32_t> tail_ubs(n_qs_in_tail_side); // lower bounds of tail
//	fmt::print("tail_size: {}\n", tail_size);
//	fmt::print("n_qs_in_tail_side: {}\n", n_qs_in_tail_side);
//	fmt::print("n_qs_in_tail_row: {}\n", n_qs_in_tail_row);
//
//	// compute the number of quadrants for each tail row by finding the max, min of that row
//#pragma omp parallel for schedule(static, 1) reduction(+:n_qs_tail)
//	for (uint32_t row_idx = wing_width; row_idx < n; row_idx += q_side_len) {
//		uint32_t max_dest_in_row = 0;
//		uint32_t min_dest_in_row = -1;
//		for (uint32_t u = row_idx; u < row_idx + q_side_len && u < n; ++u) {
//			uint64_t start = out_csr.index[u];
//			uint64_t end = out_csr.index[u + 1];
//			if (start == end || end < start) continue;
//			uint32_t max_neigh = 0;
//			uint32_t min_neigh = -1;
//			for (uint32_t offset = start; offset < end; ++offset) {
//				uint32_t v = out_csr.neighbours[offset];
//				if (v < wing_width) continue;
//				if (v >= wing_width) {
//					min_neigh = v;
//					break;
//				}
//			}
//			for (uint32_t offset = start; offset < end; ++offset) {
//				uint32_t v = out_csr.neighbours[offset];
//				if (v < wing_width) continue;
//				if (v >= max_neigh) max_neigh = v;
//			}
//			if (max_neigh > max_dest_in_row) max_dest_in_row = max_neigh;
//			if (min_neigh < min_dest_in_row) min_dest_in_row = min_neigh;
//		}
//
//
//		uint32_t tail_quad_row_length = max_dest_in_row - min_dest_in_row;
//		uint32_t n_quads_in_row = (tail_quad_row_length + q_side_len - 1) / q_side_len;
////		fmt::print("{} {} {} {}\n", min_dest_in_row, max_dest_in_row, tail_quad_row_length, n_quads_in_row);
//		uint32_t tail_offset_idx = (row_idx - wing_width) / q_side_len;
//		tail_lbs[tail_offset_idx] = min_dest_in_row;
//		tail_ubs[tail_offset_idx] = max_dest_in_row;
//		n_qs_in_tail_row[tail_offset_idx] = n_quads_in_row;
//		n_qs_tail += n_quads_in_row;
//	}
//	fmt::print("tail_lbs: {}\n", tail_lbs);
//	fmt::print("tail_ubs: {}\n", tail_ubs);
//	fmt::print("n_qs_tail: {}\n", n_qs_tail);
//	std::vector<uint32_t> cumulative_n_qs_per_row_in_tail(n_qs_in_tail_side + 1);
//	fmt::print("n_qs_in_tail_row: {}\n", n_qs_in_tail_row);
//	// now that we know the number of quadrants in the right wing, create an array of quadrants
//	std::inclusive_scan(
//		dpl::execution::par_unseq,
//		n_qs_in_tail_row.begin(),
//		n_qs_in_tail_row.end(),
//		cumulative_n_qs_per_row_in_tail.begin() + 1,
//		[](uint32_t l, uint32_t r) -> uint32_t { return l + r; },
//		0
//	);
//	fmt::print("cumulative_n_qs_per_row_in_tail: {}\n", cumulative_n_qs_per_row_in_tail);


//	Quad *tail_qs = new Quad[n_qs_tail]();
//	uint32_t *tail_q_min = new uint32_t[n_qs_tail](); // the minimum id vertex in tail quad
//	uint32_t *tail_q_max = new uint32_t[n_qs_tail](); // the maximum id vertex in tail quad
//
//	// iterate over the quadrant rows that make up the tail
//	// for each quadrant row, iterate over the range defined by the tail's lower, upper bounds
//	// and compute the number of edges contained within that quadrants to preallocate the
//	// flattened edges array
//	uint64_t total_n_edges_in_tail = 0;
//	for (uint32_t i = 0; i < n_qs_in_tail_side; ++i) {
//		uint32_t lb = tail_lbs[i];
//		uint32_t ub = tail_ubs[i];
//		uint32_t row_start = q_side_len * i + wing_width;
//		uint32_t row_end = q_side_len * (i + 1) + wing_width;
//		uint32_t quad_row = (row_start - wing_width) / q_side_len;
//		uint32_t quad_col_offset = lb - wing_width;
//		for (uint32_t col_start = lb; col_start < ub; col_start += q_side_len) {
//			uint32_t col_end = col_start + q_side_len;
////			uint32_t quad_col = (col_start - wing_width) / q_side_len;
//			uint32_t quad_col = col_start - wing_width;
//			fmt::print("quad_col_offset: {}\n", quad_col_offset);
//			fmt::print("{} {} {} {}\n", row_start, row_end, col_start, col_end);
//			fmt::print("quad_row, quad_col: {} {} \n", quad_row, quad_col);
//			uint32_t quad_idx = cumulative_n_qs_per_row_in_tail[quad_row] + \
//                          (quad_col - quad_col_offset) / q_side_len;
//			Quad &q = tail_qs[quad_idx];
//			q.q_idx = quad_idx;
//			q.qx = quad_row;
//			q.qy = quad_col;
//			uint32_t total_n_neighs_in_range = 0;
//			uint32_t min_in_quad = -1;
//			uint32_t max_in_quad = 0;
//			for (uint32_t u = row_start; u < row_end && u < n; ++u) {
//				uint32_t n_neighs_in_range = n_out_neighbours_between(u, col_start, col_end, out_csr);
//				uint32_t min_id = get_min_id_in_range(u, col_start, col_end, out_csr);
//				uint32_t max_id = get_max_id_in_range(u, col_start, col_end - 1, out_csr);
//				total_n_neighs_in_range += n_neighs_in_range;
//				if (min_id < min_in_quad) min_in_quad = min_id;
//				if (max_id > max_in_quad) max_in_quad = max_id;
//			}
//			q.nnz = total_n_neighs_in_range;
////			fmt::print("({}, {}, {})\n", min_in_quad, max_in_quad, total_n_neighs_in_range);
//			total_n_edges_in_tail += total_n_neighs_in_range;
//			tail_q_min[quad_idx] = min_in_quad;
//			tail_q_max[quad_idx] = max_in_quad;
//		}
//
//	}
//	for (int i = 0; i < n_qs_tail; ++i) {
//		fmt::print("{} {} {}\n", tail_qs[i].qx, tail_qs[i].qy, tail_qs[i].nnz);
//		fmt::print("{} {}\n", tail_q_min[i], tail_q_max[i]);
//		fmt::print("\n");
//	}
//	fmt::print("total_n_edges_in_tail: {}\n", total_n_edges_in_tail);