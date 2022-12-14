//
// Created by atrostan on 12/09/22.
//



#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>

#include "pr_experiments.h"

#include <getopt.h>
#include "typedefs.h"
#include "io.h"
#include "util.h"
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>    // std::random_shuffle
#include <cstdlib>      // std::rand, std::srand
#include <cassert>
#include <igraph.h>
#include <boost/filesystem/path.hpp>
#include "PageRank.h"
#include "sql.h"
#include "pvector.h"

int main(int argc, char *argv[]) {
	opterr = 0;
	int opt;
	ul num_vertices;
	ull num_edges;
	float percent;
	bool directed = false;
	bool debug = false;
	int num_iters = 0;
	int num_expts = 0;
	std::string input_dir;
	std::string vorder_str;
	std::string eorder_str;
	std::string sqlite_db_path;
	std::string results_path;

	while ((opt = getopt(argc, argv, "den:m:i:x:p:g:b:o:r:s:")) != -1) {
		switch (opt) {
			case 'd':
				directed = !directed;
				break;
			case 'e':
				debug = !debug;
				break;
			case 'n':
				num_vertices = atol(optarg);
				break;
			case 'm':
				num_edges = atoll(optarg);
				break;
			case 'i':
				num_iters = atoi(optarg);
				break;
			case 'x':
				num_expts = atoi(optarg);
				break;
			case 'p':
				percent = atof(optarg);
				break;
			case 'g':
				input_dir = optarg;
				break;
			case 'b':
				sqlite_db_path = optarg;
				break;
			case 'o':
				vorder_str = optarg;
				break;
			case 'r':
				results_path = optarg;
				break;
			case 's':
				eorder_str = optarg;
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
	std::string binary_graph_path = fmt::format("{}/comp.bin", input_dir);
	std::string text_graph_path = fmt::format("{}/comp.net", input_dir);

	std::string order_path = fmt::format("{}/{}", input_dir, vorder_str);
	std::string binary_order_path = fmt::format("{}/{}.bin", input_dir, vorder_str);

	std::string pr_path = fmt::format("{}/pr", input_dir);
	std::vector<double> correct_vals;
	read_binary_container<std::vector<double>>(pr_path, correct_vals);
	boost::filesystem::path p(input_dir);
	std::string graph_name = p.filename().string();

	// initializations
	std::vector<std::pair<ul, ul>> orig_edges(num_edges);
	std::vector<Edge> mapped_edges(num_edges);
	std::vector<uint32_t> iso_map;
	iso_map.resize(num_vertices, -1);

	// write header row to results csv
	std::ofstream outfile(results_path);
	outfile << "graph_name," <<
	        "datetime," <<
	        "expt_num," <<
	        "num_iters," <<
	        "vertex_order," <<
	        "edge_order," <<
	        "runtime," <<
	        "valid\n";
	outfile.close();

	// read binary edge list
	read_binary_edge_list(binary_graph_path, orig_edges);
	std::vector<Edge> indexed_edges(num_edges);

# pragma omp parallel for schedule(static)
	for (ull i = 0; i < num_edges; ++i) {
		indexed_edges[i].source = orig_edges[i].first;
		indexed_edges[i].dest = orig_edges[i].second;
		indexed_edges[i].idx = i;
	}

	std::vector<std::pair<ul, ul>>().swap(orig_edges); // delete orig edges

	if (vorder_str == "orig") {  // construct the identity map
# pragma omp parallel for schedule(static)
		for (uint32_t i = 0; i < num_vertices; ++i) { iso_map[i] = i; }
	} else {  // read (text) isomorphism map
		// read_map(order_path, iso_map);
		read_binary_container<std::vector<uint32_t>>(binary_order_path, iso_map);
	}

	iso_map.shrink_to_fit();
	assert(iso_map.size() == num_vertices);

	// translate in parallel edge list using isomorphism map
	par_translate_edge_list(indexed_edges, mapped_edges, iso_map, num_edges);

	std::string order_strings[4] = {"row", "column", "hilbert", "fgf"};
	std::map<std::string, int> edge_order_ids = {
		{"row",     0},
		{"column",  1},
		{"hilbert", 2},
		{"fgf",     3},
	};

	// sort edges according to row, column major, or hilbert  todo slashburn fgf
	const int num_edge_orders = End - Row;
	int edge_order_id = edge_order_ids[eorder_str];
	Order ord = static_cast<Order>(edge_order_id);
	fmt::print("eorder_str: {}\n", eorder_str);
	par_sort_edges(mapped_edges, ord, num_vertices);

		if (debug) {
			std::vector<std::pair<ul, ul>> sorted_edges;
			// now the edgelist is sorted, remove unneeded idx
			std::for_each(begin(mapped_edges), end(mapped_edges), [&sorted_edges](Edge e) {
				sorted_edges.emplace_back(e.source, e.dest);
			});
			// write each sorted edge list to file
			std::string eorder_path = fmt::format("{}/{}.{}", input_dir, vorder_str, eorder_str);
			write_text_edge_list(eorder_path, sorted_edges);
		}

	// run pr experiments
	for (int expt_idx = 0; expt_idx < num_expts; ++expt_idx) {

		PageRank pr = PageRank(num_vertices, num_iters, mapped_edges, percent);

		pr.compute();

		pvector<double> mapped_results(num_vertices);
#pragma omp parallel for
		for (uint32_t i = 0; i < num_vertices; ++i) {
				mapped_results[i] = pr.scores[iso_map[i]];
		}
		
		bool valid_pr = std::equal(dpl::execution::par_unseq,
		                           mapped_results.begin(), mapped_results.end(),
		                           correct_vals.begin(),
		                           [](double value1, double value2) {
			                           constexpr double epsilon = 1e-4;
			                           return std::fabs(value1 - value2) < epsilon;
		                           });
		// assert(valid_pr);
//		fmt::print("vorder_str, eorder_str: {}, {}\n", vorder_str, eorder_str);
		// fmt::print("valid_pr: {}\n", valid_pr);

		auto t = std::time(nullptr);
		auto tm = *std::localtime(&t);



		std::ostringstream oss;
		oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
		auto datetime_str = oss.str();

		PRExptRow r = {
			graph_name,
			datetime_str,
			expt_idx,
			num_iters,
			vorder_str,
			eorder_str,
			pr.runtime,
			valid_pr
		};
		write_row_to_csv(r, results_path);
		// insert_or_ignore_into_pr_expts(r, sqlite_db_path);

	}


}