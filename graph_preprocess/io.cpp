//
// Created by atrostan on 19/08/22.
//

#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include <set>
#include <fstream>
#include <map>
#include <sstream>
#include <iostream>
#include "io.h"
#include <boost/algorithm/string.hpp>
#include <chrono>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/utility.hpp>

void write_permutation(std::string path, std::map<ul, ul> &map, ul n, ull m) {
	std::ofstream outfile(path);
	outfile << fmt::format("{}\n", n);
	outfile << fmt::format("{}\n", m);

	for (auto &kv: map) {
		outfile << fmt::format("{} {}\n", kv.first, kv.second);
	}

	outfile.close();
}


void write_text_edge_list(std::string path, std::vector<std::pair<ul, ul>> &edges) {
	std::ofstream outfile(path);
// write the compressed, simplified edge list to file
	for (auto &kv: edges) {
		outfile << fmt::format("{} {}\n", kv.first, kv.second);
	}
	outfile.close();
}


void write_binary_edge_list(std::string path, std::vector<std::pair<ul, ul>> &edges) {
	std::ofstream ofs(path, std::ios::binary);
	boost::archive::binary_oarchive oa(ofs);
	oa << edges;
}

void
read_binary_edge_list_into_igraph(std::string path, std::vector<std::pair<ul, ul>> &edges, igraph_t *g, ul n, ull m,
                                  int directed) {
	std::ifstream ifs(path, std::ios::binary);
	boost::archive::binary_iarchive ia(ifs);
	ia >> edges;
	ull i = 0;
	igraph_vector_t es;
	igraph_vector_init(&es, m * 2);
	for (auto &kv: edges) {
		fmt::print("{} {}\n", kv.first, kv.second);
		VECTOR(es)[i] = kv.first;
		VECTOR(es)[i + 1] = kv.second;
		i += 2;
	}
	igraph_create(g, &es, n, 0);
	igraph_vector_destroy(&es);
}

void read_binary_edge_list(std::string path, std::vector<ul> &edges) {
	// edges' size has been preallocated and initialized
	std::ifstream ifs(path, std::ios::binary);
	boost::archive::binary_iarchive ia(ifs);
	ia >> edges;
}

void read_binary_edge_list(std::string path, std::vector<std::pair<ul, ul>> &edges) {
	// edges' size has been preallocated and initialized
	std::ifstream ifs(path, std::ios::binary);
	boost::archive::binary_iarchive ia(ifs);
	ia >> edges;
}


void read_text_edge_list(std::string path, std::vector<ul> &edges) {
	std::string line;
	ull i = 0;
	// flat edges' size has been preallocated and initialized
	std::ifstream input_file(path);
	if (input_file.is_open()) {
		while (getline(input_file, line)) {
			std::stringstream linestream(line);
			ul src;
			ul dest;
			linestream >> src >> dest;
			edges[i] = src;
			edges[i + 1] = dest;
			i += 2;
		}
		input_file.close();
	} else std::cout << "Unable to open file";
}


void read_text_edge_list(std::string path, std::vector<std::pair<ul, ul>> &edges) {
	std::string line;
	ull i = 0;
	// flat edges' size has been preallocated and initialized
	std::ifstream input_file(path);
	if (input_file.is_open()) {
		while (getline(input_file, line)) {
			std::stringstream linestream(line);
			ul src;
			ul dest;
			linestream >> src >> dest;
			edges[i].first = src;
			edges[i].second = dest;
			i += 1;
		}
		input_file.close();
	} else std::cout << "Unable to open file";
}


void write_edge_list(std::string path, std::vector<std::pair<ul, ul>> &edges, io_mode &mode) {
	std::string output_path;
	std::string file_extension;
	switch (mode) {
		case binary:
			file_extension = "bin";
			output_path = fmt::format("{}.{}", path, file_extension);
			write_binary_edge_list(output_path, edges);
			return;
		case text:
			file_extension = "net";
			output_path = fmt::format("{}.{}", path, file_extension);
			write_text_edge_list(output_path, edges);
			return;
	}
}


bool is_comment(std::string line) {
	return (line.find("%") != std::string::npos) || (line.find("#") != std::string::npos);
}

std::pair<ul, ull> read_edge_list(std::string input_path, ull m, std::vector<std::pair<ul, ul>> &edges,
                                  std::vector<std::pair<ul, ul>> &mapped_edges, std::string graph_name,
                                  std::string sqlite_db_path) {
	std::string line;
	std::ifstream input_file(input_path);
	std::map<ul, ul> iso_map;
	std::vector<ul> flat_edges(m * 2); // flattened edge list vector
//	 todo check that the number of uncommented lines in the edgelist files == m
	ull i = 0;
	ull j = 0;

	// Step 1: INGEST
	auto start = std::chrono::high_resolution_clock::now();
	if (input_file.is_open()) {
		while (getline(input_file, line)) {
			std::stringstream linestream(line);
			ul src;
			ul dest;
			linestream >> src >> dest;
			if (src == dest) continue;
			flat_edges[i] = src;
			flat_edges[i + 1] = dest;
			edges[j] = std::make_pair(src, dest);
			j += 1;
			i += 2;
		}
		input_file.close();
	} else std::cout << "Unable to open file\n";
	// truncate the flattened edge vector so that it contains only valid edges ingested
	flat_edges.resize(i);
	edges.resize(j);
	mapped_edges.resize(j);

	auto end = std::chrono::high_resolution_clock::now();
	auto ingest_time = duration_cast<time_unit>(end - start);
	// Step 2: Sort Vertices
	start = std::chrono::high_resolution_clock::now();
	// sort and remove duplicates to identify the unique vertex ids in the graph
	std::sort(dpl::execution::par_unseq, flat_edges.begin(), flat_edges.end());
	end = std::chrono::high_resolution_clock::now();
	auto sort_vs_time = duration_cast<time_unit>(end - start);
	// Step 3: Unique Vertex IDs
	start = std::chrono::high_resolution_clock::now();
	flat_edges.erase(std::unique(dpl::execution::par_unseq, flat_edges.begin(), flat_edges.end()), flat_edges.end());
	end = std::chrono::high_resolution_clock::now();
	auto unique_vs_time = duration_cast<time_unit>(end - start);

	ul max_uncompressed_vid = flat_edges.back();
	// construct a vector whose size is the max vertex id in the uncompressed repr
	// this will act as the map between uncompressed vertex ids -> vertex ids
	std::vector<ul> p(max_uncompressed_vid + 1, -1);
	ul new_id = 0;
	for (auto &orig_vid: flat_edges) {
		p[orig_vid] = new_id;
		new_id += 1;
	}
	// remap the edges
	ull eid = 0;
	for (auto &edge: edges) {
		ul src = edge.first;
		ul dest = edge.second;
		mapped_edges[eid] = std::make_pair(p[src], p[dest]);
		++eid;
	}

	// Step 4: Sort Edges
	start = std::chrono::high_resolution_clock::now();
	// sort the isomorphic edges and remove duplicates
	std::sort(dpl::execution::par_unseq, mapped_edges.begin(), mapped_edges.end());
	end = std::chrono::high_resolution_clock::now();
	auto sort_es_time = duration_cast<time_unit>(end - start);

	// Step 5: Unique Edges
	start = std::chrono::high_resolution_clock::now();
	mapped_edges.erase(std::unique(dpl::execution::par_unseq, mapped_edges.begin(), mapped_edges.end()),
	                   mapped_edges.end());
	end = std::chrono::high_resolution_clock::now();

	auto unique_es_time = duration_cast<time_unit>(end - start);

	std::vector<time_unit> times = {
			ingest_time,
			sort_vs_time,
			unique_es_time,
			sort_es_time,
			unique_es_time,
	};

	std::vector<std::string> col_labels = {
			"ingest",
			"sort_vs",
			"unique_vs",
			"sort_es",
			"unique_es",
	};


	insert_graph_into_preproc_table(graph_name, sqlite_db_path, times, col_labels);
	return std::make_pair(flat_edges.size(), mapped_edges.size());
}

void insert_graph_into_preproc_table(std::string graph_name, std::string sqlite_db_path, std::vector<time_unit> times,
                                     std::vector<std::string> col_labels) {
	sqlite3 *db;
	sqlite3_stmt *st;
	std::stringstream ss_cols;
	std::stringstream ss_vals;

	ss_cols << "(graph_name,";
	ss_vals << "(?,";

	for (int i = 0; i < col_labels.size() - 1; ++i) {
		ss_cols << col_labels[i] << ",";
		ss_vals << "?,";
	}
	ss_cols << col_labels[col_labels.size() - 1];
	ss_vals << "?)";
	ss_cols << ")";
	std::string column_label_str = ss_cols.str();
	std::string vals_str = ss_vals.str();
	std::string sql = fmt::format("INSERT INTO preproc {} VALUES {}", column_label_str, vals_str);
//	fmt::print("sql: {}\n", sql);
	if (sqlite3_open(sqlite_db_path.c_str(), &db) == SQLITE_OK) {
		sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
		sqlite3_bind_text(st, 1, graph_name.c_str(), graph_name.length(), SQLITE_TRANSIENT);

		int sql_idx = 2;
		for (auto &time: times) {
			sqlite3_bind_int(st, sql_idx, time.count());
			++sql_idx;
		}
	}
	sqlite3_step(st);
	sqlite3_finalize(st);
	sqlite3_close(db);
}


void single_val_set_int(const std::string sqlite_db_path, std::string col_name, std::string table_name,
                        std::string graph_name, int val) {

	sqlite3 *db;
	sqlite3_stmt *st;

	std::string sql = fmt::format("update {} set {} = ? where graph_name = ?", table_name, col_name);
//	fmt::print("sql: {}\n", sql);

	if (sqlite3_open(sqlite_db_path.c_str(), &db) == SQLITE_OK) {
		sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
		sqlite3_bind_int(st, 1, val);
		sqlite3_bind_text(st, 2, graph_name.c_str(), graph_name.length(), SQLITE_TRANSIENT);
	}


	sqlite3_step(st);
	sqlite3_finalize(st);
	sqlite3_close(db);
}