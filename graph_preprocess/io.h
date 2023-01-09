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
#include "igraph.h"
#include <fstream>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/utility.hpp>
#include "Quad.h"

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


void write_edge_list(std::string path, std::vector<std::pair<ul, ul>> &edges, io_mode &mode);

void write_binary_edge_list(std::string path, std::vector<std::pair<ul, ul>> &edges);

void write_text_edge_list(std::string path, std::vector<std::pair<ul, ul>> &edges);

void read_map(std::string in_path, std::vector<ul> &mp);

void write_row_to_csv(PRExptRow &r, std::string csv_path);

void read_text_degree_file(std::string path, std::vector<uint32_t> &deg);

template<typename C>
void read_binary_container(std::string inpath, C &container) {
	std::ifstream ifs(inpath, std::ios::binary);
	boost::archive::binary_iarchive ia(ifs);
	ia >> container;
}

template<typename C>
void write_binary_container(std::string outpath, C &container) {
	std::ofstream ofs(outpath, std::ios::binary);
	boost::archive::binary_oarchive oa(ofs);
	oa << container;
}

void
write_quad_array(std::string path, Quad *qs, uint32_t n_quads, uint32_t q_side_len, uint32_t wing_width, uint32_t n,
                 uint64_t m, bool is_rect, bool right_wing,
                 uint32_t n_stripes_in_right_wing,
                 uint32_t* cumulative_n_qs_per_rw_stripe);


std::vector<uint64_t> read_quad_array(std::string path, Quad *&qs, bool is_rect, bool right_wing, uint32_t *& cumulative_n_qs_per_right_wing_stripe);
namespace Eigen {
// https://scicomp.stackexchange.com/a/21438
	template<class SparseMatrix>
	inline void write_binary_sparse(const std::string &filename, const SparseMatrix &matrix) {
		assert(matrix.isCompressed() == true);
		std::ofstream out(filename, std::ios::binary | std::ios::out | std::ios::trunc);
		if (out.is_open()) {
			typename SparseMatrix::Index rows, cols, nnzs, outS, innS;
			rows = matrix.rows();
			cols = matrix.cols();
			nnzs = matrix.nonZeros();
			outS = matrix.outerSize();
			innS = matrix.innerSize();

			out.write(reinterpret_cast<char *>(&rows), sizeof(typename SparseMatrix::Index));
			out.write(reinterpret_cast<char *>(&cols), sizeof(typename SparseMatrix::Index));
			out.write(reinterpret_cast<char *>(&nnzs), sizeof(typename SparseMatrix::Index));
			out.write(reinterpret_cast<char *>(&outS), sizeof(typename SparseMatrix::Index));
			out.write(reinterpret_cast<char *>(&innS), sizeof(typename SparseMatrix::Index));

			typename SparseMatrix::Index sizeIndexS = static_cast<typename SparseMatrix::Index>(sizeof(typename SparseMatrix::StorageIndex));
			typename SparseMatrix::Index sizeScalar = static_cast<typename SparseMatrix::Index>(sizeof(typename SparseMatrix::Scalar));

//			auto p1 = matrix.outerIndexPtr();
//			for (uint32_t i = 0; i < outS; ++i) {
//				fmt::print("outerIndexPtr.: {}\n", *p1++);
//			}
////			fmt::print("matrix: {}\n", );
//			auto p2 = matrix.innerIndexPtr();
//			for (uint32_t i = 0; i < nnzs; ++i) {
//				fmt::print("innerIndexPtr.: {}\n", *p2++);
//			}

			// value is redundant - don't care about weights of edges
//			out.write(reinterpret_cast<const char*>(matrix.valuePtr()),       sizeScalar * nnzs);
			out.write(reinterpret_cast<const char *>(matrix.outerIndexPtr()), sizeIndexS * outS);
			out.write(reinterpret_cast<const char *>(matrix.innerIndexPtr()), sizeIndexS * nnzs);

			out.close();
		} else {
			std::cout << "Can not write to file: " << filename << std::endl;
		}
	}

	template<class SparseMatrix>
	inline void read_binary_sparse(const std::string &filename, SparseMatrix &matrix) {
		std::ifstream in(filename, std::ios::binary | std::ios::in);
		if (in.is_open()) {
			typename SparseMatrix::Index rows, cols, nnz, inSz, outSz;
			typename SparseMatrix::Index sizeScalar = static_cast<typename SparseMatrix::Index>(sizeof(typename SparseMatrix::Scalar));
			typename SparseMatrix::Index sizeIndex = static_cast<typename SparseMatrix::Index>(sizeof(typename SparseMatrix::Index));
			typename SparseMatrix::Index sizeIndexS = static_cast<typename SparseMatrix::Index>(sizeof(typename SparseMatrix::StorageIndex));
			std::cout << sizeScalar << " " << sizeIndex << std::endl;
			in.read(reinterpret_cast<char *>(&rows ), sizeIndex);
			in.read(reinterpret_cast<char *>(&cols ), sizeIndex);
			in.read(reinterpret_cast<char *>(&nnz  ), sizeIndex);
			in.read(reinterpret_cast<char *>(&outSz), sizeIndex);
			in.read(reinterpret_cast<char *>(&inSz ), sizeIndex);

			matrix.resize(rows, cols);
			matrix.makeCompressed();
			matrix.resizeNonZeros(nnz);

			in.read(reinterpret_cast<char *>(matrix.valuePtr()), sizeScalar * nnz);
			in.read(reinterpret_cast<char *>(matrix.outerIndexPtr()), sizeIndexS * outSz);
			in.read(reinterpret_cast<char *>(matrix.innerIndexPtr()), sizeIndexS * nnz);

			matrix.finalize();
			in.close();
		} // file is open
		else {
			std::cout << "Can not open binary sparse matrix file: " << filename << std::endl;
		}
	}
}

#endif //GRAPH_SIMPLIFY_IO_H
