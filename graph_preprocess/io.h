//
// Created by atrostan on 19/08/22.
//

#include <stdio.h>
#include <omp.h>
#include <vector>
#include "fmt/core.h"
#include <vector>
#include <set>
#include <fmt/core.h>
#include <sqlite3.h>
#include <fmt/ranges.h>
#include "typedefs.h"
#include <map>
#include <boost/serialization/vector.hpp>
#include "igraph.h"

#ifndef GRAPH_SIMPLIFY_IO_H
#define GRAPH_SIMPLIFY_IO_H


using namespace std;

using namespace fmt;


enum io_mode {
	binary,
	text
};


int segment_read(char *buff, const int len, const int count);

bool is_comment(std::string line);

int
insert_into_sqlite_db(std::string table, std::string col, ul val, std::string graph_name, std::string sqlite_db_path,
                      std::vector<time_unit> times, std::vector<std::string> col_labels);

std::pair<ul, ull>
par_read_edge_list(std::string input_path, std::vector<std::pair<ul, ul>> &mapped_edges, std::string graph_name,
                   std::string sqlite_db_path);

std::pair<ul, ull> read_edge_list(std::string input_path, ull m, std::vector<std::pair<ul, ul>> &edges,
                                  std::vector<std::pair<ul, ul>> &mapped_edges, std::string graph_name,
                                  std::string sqlite_db_path);


template<typename T>
void read_edge_list_by_mode(std::string path, std::vector<T> &edges, io_mode &mode) {
	std::string input_path;
	switch (mode) {
		case binary:
			input_path = fmt::format(path, ".bin");
			read_binary_edge_list(path, edges);
			return;
		case text:
			input_path = fmt::format(path, ".net");
			read_text_edge_list(path, edges);
			return;
	}
}


void read_binary_edge_list(std::string path, std::vector<std::pair<ul, ul>> &edges);

void read_binary_edge_list(std::string path, std::vector<ul> &edges);


void
read_binary_edge_list_into_igraph(std::string path, std::vector<std::pair<ul, ul>> &edges, igraph_t *g, ul n, ull m,
                                  int directed);

void read_text_edge_list(std::string path, std::vector<std::pair<ul, ul>> &edges);

void read_text_edge_list(std::string path, std::vector<ul> &edges);

void insert_graph_into_preproc_table(std::string graph_name, std::string sqlite_db_path, std::vector<time_unit> times,
                                     std::vector<std::string> col_labels);

void write_permutation(std::string path, std::map<ul, ul> &map, ul n, ull m);


void single_val_set_int(const std::string sqlite_db_path, std::string col_name, std::string table_name,
                        std::string graph_name, int val);

void write_edge_list(std::string path, std::vector<std::pair<ul, ul>> &edges, io_mode &mode);

void write_binary_edge_list(std::string path, std::vector<std::pair<ul, ul>> &edges);

void write_text_edge_list(std::string path, std::vector<std::pair<ul, ul>> &edges);

#endif //GRAPH_SIMPLIFY_IO_H
