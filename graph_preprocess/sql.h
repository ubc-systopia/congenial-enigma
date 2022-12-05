//
// Created by atrostan on 12/09/22.
//

#include "typedefs.h"
#include <sqlite3.h>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
//#include "fmt/core.h"

#ifndef GRAPH_PREPROCESS_SQL_H
#define GRAPH_PREPROCESS_SQL_H

void insert_graph_into_preproc_table(std::string graph_name, std::string sqlite_db_path, std::vector<time_unit> times,
                                     std::vector<std::string> col_labels);

void insert_or_ignore_into_pr_expts(PRExptRow r, std::string sqlite_db_path);

template<typename T>
void single_val_set(const std::string sqlite_db_path, std::string col_name, std::string table_name,
                    const std::string &graph_name, T val) {
	sqlite3 *db;
	sqlite3_stmt *st;
	std::stringstream ss;

	ss << "update " << table_name << " set " << col_name << " = ? where graph_name = ?" << std::endl;
	std::string sql = ss.str();
//	std::string sql = fmt::format("update {} set {} = ? where graph_name = ?", table_name, col_name);
	std::cout << "Executing:  " << sql << "\n";
	std::cout << "graph_name: " << graph_name << "\n";
	std::cout << "val: " << val << "\n";
	int sleep_ms = 1'000;
	if (sqlite3_open(sqlite_db_path.c_str(), &db) == SQLITE_OK) {
		std::cout << sqlite_db_path.c_str() << ": SQLITE_OK" << "\n";
		while (true) {
			std::cout << "Preparing statement..\n";
			int response = sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
			if (response == SQLITE_BUSY) {
				std::cout << "SQLITE_BUSY; Sleeping for " << sleep_ms << " ms;\n";
				sqlite3_sleep(sleep_ms);
			} else {

				// parse type of template param and use the appropriate sql bind fn
				if (std::is_same<T, int64_t>::value) { sqlite3_bind_int64(st, 1, val); }
				if (std::is_same<T, uint32_t>::value) { sqlite3_bind_int64(st, 1, val); }
				else if (std::is_same<T, double>::value) { sqlite3_bind_double(st, 1, val); }
				else if (std::is_same<T, int>::value) { sqlite3_bind_int(st, 1, val); }
//				else if (std::is_same<T, std::string>::value) { sqlite3_bind_text(st, 1, val.c_str()); }
				else {
					std::cout << "T: " << typeid(T).name() << "\n";
					std::cout << "INVALID TYPE..\n";
					break;
				}
				sqlite3_bind_text(st, 2, graph_name.c_str(), graph_name.length(), SQLITE_TRANSIENT);
				int step_ret = sqlite3_step(st);
				std::cout << "sql step return: " << step_ret << "\n";
				if (step_ret == SQLITE_DONE) { // successfully updated val
					std::cout << "Succesfully updated statement!\n";
					sqlite3_finalize(st);
					sqlite3_close(db);
					break;
				} else {  // retry
					std::cout << "Busy; Resetting Statement!\n";
					sqlite3_reset(st);
					sqlite3_sleep(sleep_ms);
					continue;
				}
			}
		}
	}

}

#endif //GRAPH_PREPROCESS_SQL_H
