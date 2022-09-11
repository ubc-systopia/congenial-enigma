//
// Created by atrostan on 08/09/22.
//

#include "sql.h"
#include <sqlite3.h>


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

//	std::string sql = fmt::format("INSERT INTO preproc {} VALUES {}", column_label_str, vals_str);
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
	std::stringstream ss;

	ss << "update " << table_name << " set " << col_name << " = ? where graph_name = ?" << std::endl;
	std::string sql = ss.str();
//	std::string sql = fmt::format("update {} set {} = ? where graph_name = ?", table_name, col_name);
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