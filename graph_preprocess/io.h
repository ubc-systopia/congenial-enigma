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

#ifndef GRAPH_SIMPLIFY_IO_H
#define GRAPH_SIMPLIFY_IO_H


using namespace std;

using namespace fmt;


int segment_read(char *buff, const int len, const int count);

void par_read(char *buffer, size_t size);

bool is_comment(std::string line);

int
insert_into_sqlite_db(std::string table, std::string col, ul val, std::string graph_name, std::string sqlite_db_path,
                      std::vector<time_unit> times, std::vector<std::string> col_labels);

std::pair<ul, ull> read_edge_list(std::string input_path, ull m, std::vector<std::pair<ul, ul>> &edges,
                                  std::vector<std::pair<ul, ul>> &mapped_edges, std::string graph_name,
                                  std::string sqlite_db_path);


void insert_graph_into_preproc_table(std::string graph_name, std::string sqlite_db_path, std::vector<time_unit> times,
                                     std::vector<std::string> col_labels);

void write_permutation(std::string path, std::map<ul, ul> &map, ul n, ull m);


void single_val_set_int(const std::string sqlite_db_path, std::string col_name, std::string table_name, std::string graph_name, int val);

#endif //GRAPH_SIMPLIFY_IO_H
