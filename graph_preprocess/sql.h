//
// Created by atrostan on 08/09/22.
//
#include "typedefs.h"
#include <sqlite3.h>
#include <string>
#include <vector>
#include <sstream>
#ifndef GRAPH_PREPROCESS_SQL_H
#define GRAPH_PREPROCESS_SQL_H

void insert_graph_into_preproc_table(std::string graph_name, std::string sqlite_db_path, std::vector<time_unit> times,
                                     std::vector<std::string> col_labels);
void single_val_set_int(const std::string sqlite_db_path, std::string col_name, std::string table_name,
                        std::string graph_name, int val);
#endif //GRAPH_PREPROCESS_SQL_H
