//
// Created by atrostan on 08/09/22.
//

#include <iostream>
#include "sql.h"

void insert_or_ignore_into_pr_expts(PRExptRow r, std::string sqlite_db_path) {
	std::vector<std::string> col_labels{
		"graph_name",
		"datetime",
		"expt_num",
		"num_iters",
		"vertex_order",
		"edge_order",
		"runtime",
	};

	sqlite3 *db;
	sqlite3_stmt *st;
	std::stringstream ss_cols;
	std::stringstream ss_vals;

	ss_cols << "(";
	ss_vals << "(";

	for (int i = 0; i < col_labels.size() - 1; ++i) {
		ss_cols << col_labels[i] << ",";
		ss_vals << "?,";
	}
	ss_cols << col_labels[col_labels.size() - 1];
	ss_vals << "?)";
	ss_cols << ")";
	std::string column_label_str = ss_cols.str();
	std::string vals_str = ss_vals.str();

	std::stringstream ss;

	ss << "INSERT OR IGNORE INTO pr_expts " << column_label_str << " VALUES " << vals_str << std::endl;
	std::string sql = ss.str();
	if (sqlite3_open(sqlite_db_path.c_str(), &db) == SQLITE_OK) {
		sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
		sqlite3_bind_text(st, 1, r.graph_name.c_str(), r.graph_name.length(), SQLITE_TRANSIENT);
		sqlite3_bind_text(st, 2, r.datetime.c_str(), r.datetime.length(), SQLITE_TRANSIENT);
		sqlite3_bind_int(st, 3, r.expt_num);
		sqlite3_bind_int(st, 4, r.num_iters);
		sqlite3_bind_text(st, 5, r.vertex_order.c_str(), r.vertex_order.length(), SQLITE_TRANSIENT);
		sqlite3_bind_text(st, 6, r.edge_order.c_str(), r.edge_order.length(), SQLITE_TRANSIENT);
		sqlite3_bind_int(st, 7, r.runtime);
	}
	sqlite3_step(st);
	sqlite3_finalize(st);
	sqlite3_close(db);
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

	std::stringstream ss;

	ss << "INSERT INTO preproc " << column_label_str << " VALUES " << vals_str << std::endl;
	std::string sql = ss.str();

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
                        const std::string &graph_name, int64_t val) {

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
			std::cout << sqlite_db_path.c_str() << ": SQLITE_OK" <<  "\n";
			while (true) {
			std::cout << "Preparing statement..\n";
			int response = sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
			if (response == SQLITE_BUSY) {
				std::cout << "SQLITE_BUSY; Sleeping for " <<  sleep_ms << " ms;\n";
				sqlite3_sleep(sleep_ms);
			} else {
				sqlite3_bind_int64(st, 1, val);
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