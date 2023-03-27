//
// Created by atrostan on 11/01/23.
//


#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/numeric>

#include <getopt.h>
#include <cstdint>
#include <vector>
#include <string>
#include <omp.h>
#include <likwid-marker.h>
#include <likwid.h>
#include "fmt/core.h"
#include "fmt/ranges.h"
#include "io.h"
#include "omp.h"
#include "util.h"
#include <bitset>
#include <limits.h>
#include "ips4o/ips4o.hpp"
#include "Quad.h"
#include "CSR.h"

/**
 * Preprocesses a graph that has been reordered using the Parallel-Slashburn vertex ordering
 * Preprocessing includes:
 * Splitting the graph either vertical or horizontal stripes of 4 x 4 groups of quadrants
 * (each quadrant's size = q_side_len x q_side_len)
 *
 * Groups within a stripe will be traversed either vertically or horizontally
 * Quadrants within a group will be traversed using a hilbert curve
 * Edges within a quadrant will be traversed using a hilbert curve
 *
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
	bool vertical = false;

	while ((opt = getopt(argc, argv, "bg:d:q:e:w:")) != -1) {
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
			case 'q':
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
	read_binary_container<std::vector<uint32_t >>(binary_order_path, iso_map);
	read_binary_container<std::vector<std::pair<uint32_t, uint32_t>>>(binary_edge_list_path, edge_list);

	n = iso_map.size();
	m = edge_list.size();
#pragma omp parallel for
	for (uint32_t i = 0; i < m; ++i) {
		uint32_t src = edge_list[i].first;
		uint32_t dest = edge_list[i].second;

		edge_list[i].first = iso_map[src];
		edge_list[i].second = iso_map[dest];
	}
	std::set<std::pair<uint32_t, uint32_t>> edge_set;
	if (debug) {
		for (uint32_t i = 0; i < m; ++i) {
			edge_set.insert({
				                edge_list[i].first,
				                edge_list[i].second,
			                });
		}
	}

	fmt::print("Sorting by source..\n");
	ips4o::parallel::sort(edge_list.begin(), edge_list.end());

	std::string out_deg_path = fmt::format("{}/{}", graph_dir, "out_degs");
	std::string in_deg_path = fmt::format("{}/{}", graph_dir, "in_degs");
	std::vector<uint32_t> out_degs(n);
	std::vector<uint32_t> in_degs(n);
	read_degs(out_deg_path, in_deg_path, out_degs, in_degs, n, iso_map);

	// create the in-csr and out-csr representation
	CSR in_csr = CSR(n, m, true);
	CSR out_csr = CSR(n, m, false);

	out_csr.par_populate(out_degs, edge_list);
	fmt::print("Sorting by dest..\n");
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
	uint64_t n_edges_wings = 0;
	uint64_t n_edges_tail = 0;

	// separate the edges into tail and wings; first count the number of
	// edges in each section
	// wings
#pragma omp parallel for reduction(+:n_edges_wings)
	for (uint32_t u = 0; u < n; ++u) {
		uint64_t start = out_csr.index[u];
		uint64_t end = out_csr.index[u + 1];
		if (u < wing_width) {
			n_edges_wings += end - start;
		} else {
			for (uint64_t i = start; i < end; ++i) {
				uint32_t v = out_csr.neighbours[i];
				if (v >= wing_width) { break; }
				++n_edges_wings;
			}
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
			++n_edges_tail;
		}
	}
	assert(n_edges_wings + n_edges_tail == m);
	fmt::print("m: {}\n", m);
	fmt::print("n: {}\n", n);
	fmt::print("n_edges_lwing: {}\n", n_edges_wings);
	fmt::print("n_edges_tail: {}\n", n_edges_tail);
	fmt::print("q_side_len: {}\n", q_side_len);
	fmt::print("wing_width: {}\n", wing_width);

	uint32_t stripe_len = q_side_len * 4;

	CSR *csr;
	vertical ? csr = &in_csr : csr = &out_csr;

	// compute the number of stripes in each wing
	uint32_t n_stripes = (n + stripe_len) / stripe_len;
	uint32_t n_qs_per_side = (n + q_side_len) / q_side_len;
	uint32_t n_qs_wings = 0;
// a slice can be either a row or col
	std::vector<uint32_t> n_qs_in_slice(n_qs_per_side);

#pragma omp parallel for schedule(static, 1) reduction(+:n_qs_wings)
	for (uint32_t slice_idx = 0; slice_idx < n; slice_idx += q_side_len) {
		uint32_t max_dest_in_row = 0;
		for (uint32_t i = slice_idx; i < slice_idx + q_side_len && i < n; ++i) {
			uint64_t start = csr->index[i];
			uint64_t end = csr->index[i + 1];
			uint32_t largest_dest_seen = 0;
			if (i >= wing_width) {

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
				largest_dest_seen = out_csr.neighbours[offset - 1];

			} else {
				largest_dest_seen = out_csr.neighbours[end - 1];
			}
			if (largest_dest_seen > max_dest_in_row) max_dest_in_row = largest_dest_seen;
		}
		uint32_t n_quads_in_slice = (max_dest_in_row + q_side_len) / q_side_len;
		n_qs_in_slice[slice_idx / q_side_len] = n_quads_in_slice;
		n_qs_wings += n_quads_in_slice;
	}

	fmt::print("n_qs_in_stripe: {}\n", n_qs_in_slice);
	std::vector<uint32_t> cumulative_n_qs_per_slice_in_wings(n_qs_per_side + 1);
	// now that we know the number of quadrants in the right wing create an array of quadrants
	std::inclusive_scan(
		dpl::execution::par_unseq,
		n_qs_in_slice.begin(),
		n_qs_in_slice.end(),
		cumulative_n_qs_per_slice_in_wings.begin() + 1,
		[](uint32_t l, uint32_t r) -> uint32_t { return l + r; },
		0
	);
	fmt::print("cumulative_n_qs_per_row_in_wings: {}\n", cumulative_n_qs_per_slice_in_wings);
	Quad *wing_qs = new Quad[n_qs_wings]();

	// assign the quad coordinates for the quadrants in the left wing
#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < n_qs_per_side; i++) {
		uint32_t n_qs = n_qs_in_slice[i];
		for (uint32_t j = 0; j < n_qs; ++j) {
			uint32_t quad_idx = cumulative_n_qs_per_slice_in_wings[i] + j;
			Quad &q = wing_qs[quad_idx];
			q.q_idx = quad_idx;
			if (vertical) {
				q.qx = j;
				q.qy = i;
			} else {
				q.qx = i;
				q.qy = j;
			}
		}
	}

	// count the number of edges that will reside in each quadrant of the wings
	// this can be done in parallel since each thread exclusively inspects a quadrant
	// slice (row or col)
#pragma omp parallel for schedule(static)
	for (uint32_t slice_idx = 0; slice_idx < n; slice_idx += q_side_len) {
		for (uint32_t i = slice_idx; i < slice_idx + q_side_len && i < n; ++i) {
			uint64_t start = csr->index[i];
			uint64_t end = csr->index[i + 1];
			uint32_t quad_i = i / q_side_len;

			if (i >= wing_width) {
				// skip vertices that do not have any neighbours that reside in the wing
				if (csr->neighbours[start] >= wing_width || start == end) { continue; }
				for (uint64_t offset = start; offset < end; offset++) {
					uint32_t j = csr->neighbours[offset];
					if (j >= wing_width) break;
					uint32_t quad_j = j / q_side_len;
					uint32_t quad_idx = cumulative_n_qs_per_slice_in_wings[quad_i] + quad_j;
					Quad &q = wing_qs[quad_idx];
					q.nnz++;
				}
			} else {
				if (start == end) { continue; }
				for (uint64_t offset = start; offset < end; offset++) {
					uint32_t j = csr->neighbours[offset];
					uint32_t quad_j = j / q_side_len;
					uint32_t quad_idx = cumulative_n_qs_per_slice_in_wings[quad_i] + quad_j;
					Quad &q = wing_qs[quad_idx];
					q.nnz++;
				}
			}
		}
	}

#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < n_qs_wings; ++i) {
		Quad &q = wing_qs[i];
		q.edges = new uint32_t[q.nnz * 2]();
	}

	uint32_t *n_edges_seen_in_wing_q = new uint32_t[n_qs_wings]();

	for (uint32_t stripe_idx = 0; stripe_idx < n_stripes; stripe_idx++) {
		uint32_t stripe_start = stripe_idx * stripe_len;
		uint32_t stripe_end = (stripe_idx + 1) * stripe_len;
#pragma omp parallel for schedule(static, 1)
		for (uint32_t row_idx = stripe_start; row_idx < stripe_end; row_idx += q_side_len) {
			for (uint32_t i = row_idx; i < row_idx + q_side_len && i < n; ++i) {
				uint64_t start = csr->index[i];
				uint64_t end = csr->index[i + 1];
				uint32_t quad_i = i / q_side_len;
				uint32_t u, v;
				if (i >= wing_width) {
					// skip vertices that do not have any neighbours that reside in the wing
					if (csr->neighbours[start] >= wing_width || start == end) { continue; }
					for (uint64_t offset = start; offset < end; offset++) {
						uint32_t j = csr->neighbours[offset];
						if (j >= wing_width) break;
						uint32_t quad_j = j / q_side_len;
						uint32_t quad_idx = cumulative_n_qs_per_slice_in_wings[quad_i] + quad_j;
						Quad &q = wing_qs[quad_idx];
						vertical ? u = j, v = i : u = i, v = j;
						auto p = edge_offset(u, v, q_side_len);
						uint32_t qu = p.first;
						uint32_t qv = p.second;
						// insert the offset (u, v) to the flattened edges array of q
						uint32_t edge_offset_into_q = n_edges_seen_in_wing_q[quad_idx];
						q.edges[edge_offset_into_q] = qu;
						q.edges[edge_offset_into_q + 1] = qv;
						n_edges_seen_in_wing_q[quad_idx] += 2;
					}
				} else {
					if (start == end) { continue; }
					for (uint64_t offset = start; offset < end; offset++) {
						uint32_t j = csr->neighbours[offset];
						uint32_t quad_j = j / q_side_len;
						uint32_t quad_idx = cumulative_n_qs_per_slice_in_wings[quad_i] + quad_j;
						Quad &q = wing_qs[quad_idx];
						vertical ? u = j, v = i : u = i, v = j;
						auto p = edge_offset(u, v, q_side_len);
						uint32_t qu = p.first;
						uint32_t qv = p.second;
// insert the offset (u, v) to the flattened edges array of q
						uint32_t edge_offset_into_q = n_edges_seen_in_wing_q[quad_idx];
						q.edges[edge_offset_into_q] = qu;
						q.edges[edge_offset_into_q + 1] = qv;
						n_edges_seen_in_wing_q[quad_idx] += 2;
					}
				}
			}
		}
	}

	uint64_t total_edges_seen = 0;
	for (uint32_t i = 0; i < n_qs_wings; ++i) {
		Quad &q = wing_qs[i];
//		fmt::print("q.qx, q.qy: {} {}\n", q.qx, q.qy);
		total_edges_seen += q.nnz;
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
//	uint32_t n_rects_in_tail = (tail_size + q_side_len - 1) / q_side_len;
	uint32_t n_rects_in_tail = (tail_size + q_side_len) / q_side_len;
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
	uint64_t w_count_test = 0;
	uint64_t tail_count_test = 0;

#pragma omp parallel for schedule(static) reduction(+:w_count_test)
	for (uint32_t i = 0; i < n_qs_wings; ++i)
		w_count_test += wing_qs[i].nnz;

#pragma omp parallel for schedule(static) reduction(+:tail_count_test)
	for (uint32_t i = 0; i < n_rects_in_tail; ++i)
		tail_count_test += tail_rects[i].nnz;

	assert(w_count_test == n_edges_wings);
	assert(tail_count_test == n_edges_tail);
	assert(w_count_test + tail_count_test == m);

	if (debug) {
		for (uint64_t i = 0; i < n_qs_wings; ++i) {
			Quad &q = wing_qs[i];
			for (uint32_t j = 0; j < q.nnz; ++j) {
				uint32_t src = q.edges[j * 2];
				uint32_t dest = q.edges[j * 2 + 1];
				// wings
				uint32_t u = q_side_len * q.qx + src;
				uint32_t v = q_side_len * q.qy + dest;
				if (not edge_set.count({u, v})) {
					fmt::print("u, v: {}\n", u, v);
				}
			}
		}
	}


//	hilbert order the quadrants in the wings
// quadrants are chunked into rows of size stripe_len
	uint32_t n_quads_per_stripe = stripe_len / q_side_len; // both powers of 2 so even division

	uint32_t *cumulative_n_quads_per_wing_stripe = new uint32_t[n_stripes + 1]();
	uint32_t wing_stripe_idx = 1;
	for (uint32_t i = n_quads_per_stripe; i < n_qs_per_side; i += n_quads_per_stripe) {
		cumulative_n_quads_per_wing_stripe[wing_stripe_idx] = cumulative_n_qs_per_slice_in_wings[i];
		++wing_stripe_idx;
	}
	cumulative_n_quads_per_wing_stripe[n_stripes] = cumulative_n_qs_per_slice_in_wings[n_qs_per_side];

	print_arr<uint32_t>(cumulative_n_quads_per_wing_stripe, n_stripes + 1);
// hilbert sort all quadrants within a stripe
	for (uint32_t i = 0; i < n_stripes; ++i) {
		uint32_t stripe_start = cumulative_n_quads_per_wing_stripe[i];
		uint32_t stripe_end = cumulative_n_quads_per_wing_stripe[i + 1];
	}
}

