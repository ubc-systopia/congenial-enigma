#include <iostream>
#include <getopt.h>
#include <string>
#include "io.h"
#include "typedefs.h"
#include <set>
#include <algorithm>    // std::random_shuffle
#include <cstdlib>      // std::rand, std::srand
#include <boost/filesystem.hpp>
#include <eigen3/Eigen/Sparse>
#include "pvector.h"
#include <eigen3/Eigen/Core>
#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include "spray.hpp"

void convert_to_sparse_mat_and_save(bool laplacian, bool in_degree,
                                    std::vector<std::pair<uint32_t, uint32_t>> &edges,
                                    uint32_t n, std::string out_path) {
	typedef Eigen::Triplet<ul> T;
	pvector<uint32_t> deg(n, 0);
	uint64_t m = edges.size();
	uint64_t n_entries_in_mat = laplacian ? n + m : m;
	pvector<T> triples(n_entries_in_mat);
	uint32_t mat_val = laplacian ? -1 : 1;

#pragma omp parallel for schedule(static)
	for (uint64_t i = 0; i < m; ++i) {
		uint32_t src = edges[i].first;
		uint32_t dest = edges[i].second;
		triples[i] = T(src, dest, mat_val);
		if (in_degree) { deg[dest]++; }
		else { deg[src]++; }
	}
	if (laplacian) {
#pragma omp parallel for schedule(static)
		for (uint32_t i = 0; i < deg.size(); ++i) {
			triples[i + m] = T(i, i, deg[i]);
		}
	}

	Eigen::SparseMatrix<ul> M(n, n);
	M.setFromTriplets(triples.begin(), triples.end());
	write_binary_sparse(out_path, M);
}

template<typename T>
void sort_and_remove_dups(std::vector<T> &v) {
	std::sort(dpl::execution::par_unseq, v.begin(), v.end());
	v.erase(std::unique(dpl::execution::par_unseq, v.begin(), v.end()), v.end());
}

bool is_line_comment(std::string line) {
	return (line[0] == '#' || line[0] == '%');
}

typedef std::tuple<uint32_t, uint32_t, double> WeightedEdge;

void write_bipartite_to_disk(std::vector<WeightedEdge> &es, std::string outpath) {
	std::ofstream outfile(outpath);
	for (WeightedEdge &e: es) {
		outfile << fmt::format("{} {} {}\n", get<0>(e), get<1>(e), get<2>(e));
	}
	outfile.close();

}

void compress_bipartite(std::string input_path, uint64_t m, std::string output_path) {
	std::vector<uint32_t> srcs(m);
	std::vector<uint32_t> dests(m);
	std::vector<WeightedEdge> es(m);
	// read the edgelist
	std::string line;
	fmt::print("input_path: {}\n", input_path);
	std::ifstream input_file(input_path);
	// ingest
	uint64_t i = 0;
	if (input_file.is_open()) {
		while (getline(input_file, line)) {
			std::stringstream linestream(line);
			if (is_comment(line)) continue;
			uint32_t v1;
			uint32_t v2;
			double wt;
			linestream >> v1 >> v2 >> wt;
//			fmt::print("v1, v2, wt: {} {}  {}\n", v1, v2, wt);
			es[i] = {v1, v2, wt};
			srcs[i] = v1;
			dests[i] = v2;
			++i;
		}
		input_file.close();
	} else std::cout << "Unable to open file\n";
	m = i;
	fmt::print("m: {}\n", m);
	fmt::print("srcs, dests: {} {} \n", srcs.size(), dests.size());
	sort_and_remove_dups<uint32_t>(srcs);
	sort_and_remove_dups<uint32_t>(dests);
	fmt::print("srcs, dests: {} {} \n", srcs.size(), dests.size());

	// create the vertex id map to the compressed space
	uint32_t n1 = srcs.size();
	uint32_t n2 = dests.size();
	uint32_t mx_v1 = srcs[n1 - 1];
	uint32_t mv_v2 = dests[n2 - 1];
	fmt::print("mx_v1, mv_v2: {} {}\n", mx_v1, mv_v2);
	std::vector<uint32_t> V1(mx_v1);
	std::vector<uint32_t> V2(mv_v2);
	fmt::print("srcs: {}\n", srcs);
#pragma omp parallel for
	for (uint32_t i = 0; i < n1; ++i) { V1[srcs[i]] = i; }
	fmt::print("n1, n2, m: {} {} {}\n", n1, n2, m);
	// offset the vertex IDs in V2 so that their ID does not overlap with the IDs of V1
#pragma omp parallel for
	for (uint32_t i = 0; i < n2; ++i) {
		
		V2[dests[i]] = n1 + i;
	}

	fmt::print("V1: {}\n", V1);
	fmt::print("V2: {}\n", V2);

#pragma omp parallel for
	for (uint64_t i = 0; i < m; ++i) {
		WeightedEdge &e = es[i];
		get<0>(e) = V1[get<0>(e)];
		get<1>(e) = V2[get<1>(e)];
	}
	auto i1 = std::adjacent_find(es.begin(), es.end());
	uint32_t idx = std::distance(es.begin(), i1);
	fmt::print("es.size(): {}\n", es.size());
	sort_and_remove_dups<WeightedEdge>(es);
	fmt::print("es.size(): {}\n", es.size());

	//write edges to disk
	write_bipartite_to_disk(es, output_path);

	// compute the degrees of the original vertex set
	uint32_t *v1_degs = new uint32_t[n1]();
	uint32_t *v2_degs = new uint32_t[n2]();
	spray::BlockReduction4096<uint32_t> v1_p(n1, v1_degs);
	spray::BlockReduction4096<uint32_t> v2_p(n2, v2_degs);


#pragma omp parallel for reduction(+:v1_p) reduction(+:v2_p)
	for (uint64_t i = 0; i < m; ++i) {
		WeightedEdge &e = es[i];
		v1_p[get<0>(e)]++;
		v2_p[get<1>(e) - n1]++;
	}

	uint32_t *res = std::max_element(v1_degs, v1_degs + n1);
	fmt::print("std::distance(v1_degs, res): {}\n", std::distance(v1_degs, res));
	fmt::print("*res: {}\n", *res);


	uint32_t n = n1 + n2;
	std::vector<uint32_t> iso_map(n);

#pragma omp parallel for
	for (uint32_t i = 0; i < n; ++i) iso_map[i] = i;

	// sort V1, V2 by ascending degree (separately)
	std::sort(dpl::execution::par_unseq,
	          iso_map.begin(), iso_map.begin() + n1,
	          [&](uint32_t v11, uint32_t v12) -> bool {
		          return v1_degs[v11] > v1_degs[v12];
	          });
//
//	std::sort(dpl::execution::par_unseq,
//	          iso_map.begin() + n1, iso_map.begin() + n1 + n2,
//	          [&](uint32_t v21, uint32_t v22) -> bool {
//		          return v2_degs[v21 - n1] > v2_degs[v22 - n1];
//	          });

	// write the isomoprhism map

	for (uint32_t i = 0; i < 10; ++i) {
		fmt::print("iso_map[i]: {}\n", iso_map[i]);
	}
//	fmt::print("iso_map: {}\n", iso_map);
	write_isomap(fmt::format("{}.{}", output_path, "srt"), iso_map, n, n1, n2, m);
	delete[]v1_degs;
	delete[]v2_degs;
	return;
}

/**
 * Read an undirected/directed graph from an edge-list file
 * Simplify the graph (Removes loop and/or multiple edges from the graph)
 * Compress the graph's ID space
 * @return
 */

int main(int argc, char *argv[]) {
	opterr = 0;
	int opt;
	ull num_edges = 0;
	std::string input_path;
	std::string sqlite_db_path;
	std::string output_path;
	std::vector<io_mode> io_modes;
	bool bipartite = false;
	uint32_t n1 = 0;
	uint32_t n2 = 0;

	while ((opt = getopt(argc, argv, "ptig:b:m:o:u:v:")) != -1) {
		switch (opt) {
			case 'p':
				bipartite = !bipartite;
				break;
			case 't':
				io_modes.push_back(text);
				break;
			case 'i':
				io_modes.push_back(binary);
				break;
			case 'g':
				input_path = optarg;
				break;
			case 'b':
				sqlite_db_path = optarg;
				break;
			case 'm':
				num_edges = atoll(optarg);
				break;
			case 'u':
				n1 = atol(optarg);
				break;
			case 'v':
				n2 = atol(optarg);
				break;
			case 'o':
				output_path = optarg;
				break;
			case '?':
				if (optopt == 'k')
					printf("Option -%c requires a long long.\n", optopt);
				else if (optopt == 'g' || optopt == 'o')
					printf("Option -%c requires a string.\n", optopt);
				else
					printf("Unknown option character '\\x%x'.\n", optopt);
				return 1;
			default:
				abort();
		}
	}

	if (bipartite) {
		// compress a weighted bipartite graph
		// given a symmetric, weighted bipartite (e.g. user ratings), with vertex sets V1, V2, with n1, n2 vertices per
		// set, respectively, compress the vertex ID space of the vertex sets so that their vertex IDs lie within the range
		// [0, n1), [0, n2).
		compress_bipartite(input_path, num_edges, output_path);


	} else {
		// if no io mode supplied, default to text
		if (io_modes.empty()) {
			io_modes.push_back(text);
		}

		std::vector<std::pair<ul, ul>> mapped_edges;
		mapped_edges.resize(num_edges);

		boost::filesystem::path p(input_path);
		boost::filesystem::path dir = p.parent_path();
		std::string graph_name = dir.filename().string();

		std::vector<std::pair<ul, ul>> edges;
		edges.resize(num_edges);
		std::pair<ul, ull> nm = read_edge_list(input_path, num_edges, edges, mapped_edges, graph_name, sqlite_db_path);



//	std::pair<ul, ull> nm = par_read_edge_list(input_path, mapped_edges, graph_name, sqlite_db_path);
		ul n = nm.first;
		ull m = nm.second;
		fmt::print("n: {}\n", n);
		fmt::print("m: {}\n", m);

		for (auto &io_mode: io_modes) {
			write_edge_list(output_path, mapped_edges, io_mode);
		}
	}


	return 0;

}
