//
// Created by atrostan on 18/11/22.
//

#include <omp.h>
#include <map>
#include "command_line.h"
#include "cc.h"
#include "fmt/core.h"
#include "fmt/ranges.h"
#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/config.hpp>
#include "io.h"
#include <boost/filesystem.hpp>

bool in_gcc(NodeID vid, NodeID gcc_id, pvector<NodeID> &comp) {
	return comp[vid] == gcc_id;
}
void write_cc(pvector<NodeID> &comp, NodeID gcc_id, Graph &g, std::string out_path) {
	// a mapping between the original vertex to the vertex ids of vertices
	// that are a part of the Giant Connected Component
	std::unordered_map<uint32_t, uint32_t> p;
	uint32_t gcc_count = 0;
	for (uint32_t i = 0; i < g.num_nodes(); ++i) {
		if (in_gcc(i, gcc_id, comp)) {
			p[i] = gcc_count;
			++gcc_count;
		}
	}

	uint64_t n_edges_in_gcc = 0;
#pragma omp parallel for schedule(dynamic, 16384) reduction(+:n_edges_in_gcc)
	for (uint64_t u = 0; u < g.num_nodes(); u++) {
		if (!in_gcc(u, gcc_id, comp)) continue;
		for (uint64_t v: g.out_neigh(u)) {
			if (!in_gcc(v, gcc_id, comp)) continue;
			++n_edges_in_gcc;
		}
	}
	pvector<std::pair<uint32_t, uint32_t>> mapped_edges(n_edges_in_gcc);
	n_edges_in_gcc = 0;
	for (uint64_t u = 0; u < g.num_nodes(); u++) {

		if (!in_gcc(u, gcc_id, comp)) continue;
		for (uint64_t v: g.out_neigh(u)) {
			if (!in_gcc(v, gcc_id, comp)) continue;
			mapped_edges[n_edges_in_gcc] = {p[u], p[v]};
			++n_edges_in_gcc;
		}
	}
	std::sort(dpl::execution::par_unseq, mapped_edges.begin(), mapped_edges.end());
	// write sorted, mapped edges to outfile
	std::ofstream outfile(out_path);
	for (const auto &kv: mapped_edges) {
		outfile << fmt::format("{} {}\n", kv.first, kv.second);
	}
	outfile.close();
}

std::pair<NodeID, uint32_t> compute_gcc_size(pvector<NodeID> &comp, std::vector<NodeID> &unique_cids) {
	int num_threads = omp_get_max_threads();
	// thread local map of number of vertices in each component
	std::vector<std::map<NodeID, uint32_t>> counts(num_threads);
	for (uint i = 0; i < num_threads; ++i) {
		for (const auto &cid: comp) {
			counts[i][cid] = 0;
		}
	}
#pragma omp parallel for schedule(static)
	for (uint32_t i = 0; i < comp.size(); ++i) {
		++counts[omp_get_thread_num()][comp[i]];
	}

	std::map<NodeID, uint32_t> global_counts;
	for (const auto &cid: unique_cids) {
		global_counts[cid] = 0;
	}
	for (uint i = 0; i < num_threads; ++i) {
		for (const auto &kv: counts[i]) {
			global_counts[kv.first] += kv.second;
		}
	}
	std::map<NodeID, uint32_t>::iterator gcc
		= std::max_element(global_counts.begin(), global_counts.end(),
		                   [](const std::pair<NodeID, uint32_t> &a, const std::pair<NodeID, uint32_t> &b) -> bool {
			                   return a.second < b.second;
		                   });

	NodeID gcc_id = gcc->first;
	uint32_t gcc_size = gcc->second;
	return std::make_pair(gcc_id, gcc_size);
}

void scc(std::string graph_path, Graph &csr, std::string out_path) {
	using namespace boost;

	typedef adjacency_list<vecS, vecS, directedS> BoostGraph;
	std::vector<std::pair<uint32_t, uint32_t>> edges(csr.num_edges());
	read_text_edge_list(graph_path, edges);
	BoostGraph G(csr.num_nodes());
	for (auto &kv: edges) { add_edge(kv.first, kv.second, G); }
	typedef graph_traits<adjacency_list<vecS, vecS, directedS> >::vertex_descriptor Vertex;
	typedef graph_traits<adjacency_list<vecS, vecS, directedS> >::vertex_descriptor Vertex;

	pvector<NodeID> comp(num_vertices(G));
	pvector<uint32_t> discover_time(num_vertices(G));
	std::vector<default_color_type> color(num_vertices(G));
	std::vector<Vertex> root(num_vertices(G));
	int num = strong_components(G, make_iterator_property_map(comp.begin(), get(vertex_index, G)),
	                            root_map(make_iterator_property_map(root.begin(), get(vertex_index, G))).
		                            color_map(make_iterator_property_map(color.begin(), get(vertex_index, G))).
		                            discover_time_map(
		                            make_iterator_property_map(discover_time.begin(), get(vertex_index, G))));

//	std::cout <<  num << std::endl;
	std::vector<int>::size_type i;

	std::vector<NodeID> unique_cids(comp.begin(), comp.end());

	std::sort(dpl::execution::par_unseq, unique_cids.begin(), unique_cids.end());
	// get unique component ids
	auto last = std::unique(dpl::execution::par_unseq, unique_cids.begin(), unique_cids.end());
	unique_cids.erase(last, unique_cids.end());
	auto res = compute_gcc_size(comp, unique_cids);
	NodeID gcc_id = res.first;
	uint32_t gcc_size = res.second;
//	fmt::print("gcc_id, gcc_size {} {}\n", gcc_id, gcc_size);
//	fmt::print(": {}\n", );
	write_cc(comp, gcc_id, csr, out_path);
	fmt::print("{} {}\n", num, gcc_size);

	return;
}

void write_symm_degs(std::string input_path, Graph &g) {
	boost::filesystem::path p(input_path);
	boost::filesystem::path dir = p.parent_path();
	std::string graph_name = dir.filename().string();
	std::string degs_path = fmt::format("{}/degs", dir.string());
	std::ofstream degs_file(degs_path);
	for (uint32_t u = 0; u < g.num_nodes(); ++u) {
		degs_file << g.out_degree(u) << "\n";
	}
	degs_file.close();
}

void write_degs(std::string input_path, Graph &g) {
	boost::filesystem::path p(input_path);
	boost::filesystem::path dir = p.parent_path();
	std::string graph_name = dir.filename().string();
	std::string out_degs_path = fmt::format("{}/out_degs", dir.string());
	std::string in_degs_path = fmt::format("{}/in_degs", dir.string());
	std::ofstream out_degs_file(out_degs_path);
	std::ofstream in_degs_file(in_degs_path);
	for (uint32_t u = 0; u < g.num_nodes(); ++u) {
		out_degs_file << g.out_degree(u) << "\n";
		in_degs_file << g.in_degree(u) << "\n";
	}
	out_degs_file.close();
	in_degs_file.close();
}

int main(int argc, char *argv[]) {
	CLApp cli(argc, argv, "computes the largest strongly or weakly connected component of an input graph");
	if (!cli.ParseArgs())
		return -1;

	Builder b(cli);
	Graph g = b.MakeGraph();
	bool directed = !cli.symmetrize();
	if (directed) {
//		 use boost to compute strongly connected components;
		write_degs(cli.filename(), g);
		scc(cli.filename(), g, cli.out_filename());
		return 0;
	}

	write_symm_degs(cli.filename(), g);
	pvector<NodeID> comp = Afforest(g, 2);

	std::vector<NodeID> unique_cids(comp.begin(), comp.end());

	std::sort(dpl::execution::par_unseq, unique_cids.begin(), unique_cids.end());
	// get unique component ids
	auto last = std::unique(dpl::execution::par_unseq, unique_cids.begin(), unique_cids.end());
	unique_cids.erase(last, unique_cids.end());
	// todo write number of ccs, and size of gcc to db
	auto res = compute_gcc_size(comp, unique_cids);
	NodeID gcc_id = res.first;
	uint32_t gcc_size = res.second;
	fmt::print("{} {}\n", unique_cids.size(), gcc_size);

	write_cc(comp, gcc_id, g, cli.out_filename());
}