//
// Created by atrostan on 29/11/22.
//

#include <map>
#include <string>
#include <vector>
#include "fmt/core.h"
#include "fmt/ranges.h"
#include "io.h"
#include "omp.h"
#include "sql.h"
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/StdVector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>

#include <Spectra/SymEigsSolver.h>
#include "pvector.h"
#include <type_traits>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/SymGEigsShiftSolver.h>
#include <Spectra/GenEigsRealShiftSolver.h>

#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseGenRealShiftSolve.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>

#include <Spectra/contrib/PartialSVDSolver.h>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>

#ifndef GRAPH_PREPROCESS_STATS_H
#define GRAPH_PREPROCESS_STATS_H


/**
 *
 * @tparam T type of matrix entries
 * @tparam SolveOP type of eigen solve operation
 * @tparam Solver type of eigen solver
 * @tparam MatClass type of sparse matrix
 */
template<
	class T,
	template<class> class SolveOP,
	template<class> class Solver>
void solve_shift_invert(Eigen::SparseMatrix<T> &mat, Spectra::SortRule eig_sort_rule, double sigma) {
	SolveOP<T> op(mat);
	Solver<SolveOP<T>> eigs(op, 6, 10, sigma);
	eigs.init();
	eigs.compute(eig_sort_rule);
	if (eigs.info() == Spectra::CompInfo::Successful) {
		auto evals = eigs.eigenvalues();
		for (const auto &eval: evals) {
			fmt::print("eval: {}\n", eval);
			auto converted = 1 / (eval - sigma);
			fmt::print("converted: {}\n", converted);
		}

		// Will get (3.0, 2.0, 1.0)
	} else {
		fmt::print("eigs.info: {}\n", int(eigs.info()));
	}
}

template<typename T>
double solve_directed(Eigen::SparseMatrix<T> &A, int n_eig_vals) {
	Eigen::initParallel();
	Spectra::SparseGenMatProd<T> op(A);

	Spectra::GenEigsSolver<Spectra::SparseGenMatProd<T>> eigs(op, n_eig_vals, n_eig_vals * 10);
	eigs.init();
	eigs.compute(Spectra::SortRule::LargestMagn);
	if (eigs.info() == Spectra::CompInfo::Successful) {
		auto evals = eigs.eigenvalues();
		return evals[0].real();
	} else {
		return -1;
	}
}

template<typename T>
std::pair<double, double> solve_symmetric(Eigen::SparseMatrix<T> &A, int n_eig_vals) {
	Eigen::initParallel();

	Spectra::SparseSymMatProd<T> op(A);
	Spectra::SymEigsSolver<Spectra::SparseSymMatProd<T>> eigs(op, n_eig_vals, n_eig_vals * 10);
	eigs.init();
	eigs.compute(Spectra::SortRule::LargestMagn);
	if (eigs.info() == Spectra::CompInfo::Successful) {
		auto evals = eigs.eigenvalues();
		return {evals[1], evals[1] / evals[0]};
	} else {
		return {-1, -1};
	}
}

template<typename T>
double compute_svd(Eigen::SparseMatrix<T> &A, int n_vals) {
	Eigen::initParallel();

	Spectra::PartialSVDSolver<Eigen::SparseMatrix<T>> svds(A, n_vals, n_vals * 3);
	int nconv = svds.compute();
	auto svals = svds.singular_values();
	return svals[0];
}

template<typename T, template<class> class C>
void populate_mat_from_edgelist(Eigen::SparseMatrix<T> &A, std::string in_path, uint32_t n, bool laplacian, T mat_val) {
	fmt::print("in_path: {}\n", in_path);
	C<std::pair<uint32_t, uint32_t>> edges;
	read_binary_container(in_path, edges);
	fmt::print("read_binary_container.\n");
	typedef Eigen::Triplet<T> Triplet;
	pvector<uint32_t> deg(n, 0);
	uint64_t m = edges.size();
	uint64_t n_entries_in_mat = laplacian ? n + m : m;
	pvector<Triplet> triples(n_entries_in_mat);

#pragma omp parallel for schedule(static)
	for (uint64_t i = 0; i < m; ++i) {
		uint32_t src = edges[i].first;
		uint32_t dest = edges[i].second;
		triples[i] = Triplet(src, dest, mat_val);
	}
	fmt::print("populated triples with edges.\n");

	if (laplacian) {
		// add 1 to the diagonal (to store the degree of the vertex later)
#pragma omp parallel for schedule(static)
		for (uint32_t i = 0; i < n; ++i) {
			triples[i + m] = Triplet(i, i, 1);
		}
	}

	fmt::print("Creating Sparse Mat from triples..\n");
	A.setFromTriplets(triples.begin(), triples.end());

}

/**
 *
 * @tparam T type of matrix entries; int32_t for laplacian, uint32_t otherwise
 * @tparam C type of edge list container - either pvector, or std::vector of std::pair<uint32_t, uint32_t>
 * @param symmetric Whether the matrix should be symmetrized
 * @param laplacian Should we compute the laplacian
 * @param in_path Path to the binary file of the graph's edgelist
 * @param n number of nodes in the graph
 * @param matval value to store at A_ij to indicate the existence of the edge i -> j
 */
template<class T, template<class> class C>
void compute_eig_stats(bool symmetric, bool laplacian, std::string in_path, uint32_t n, T mat_val,
                       std::string sqlite3_db_path) {


	Eigen::SparseMatrix<T> A(n, n);
	populate_mat_from_edgelist<T, C>(A, in_path, n, laplacian, mat_val);
	fmt::print("Created.\n");
	fmt::print("computing reciprocity..\n");
	uint64_t recip_edges = 0;
#pragma omp parallel for schedule(dynamic) reduction(+:recip_edges)
	for (int k = 0; k < A.outerSize(); ++k)
		for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, k); it; ++it) {
			// if the opposite edge exists
			if (A.coeff(k, it.index()) == 1) {
				++recip_edges;
			}
		}

	// save adj list as binary for future numpy read
	boost::filesystem::path pth(in_path);
	boost::filesystem::path dir = pth.parent_path();
	std::string graph_name = dir.filename().string();
	std::string mat_path = fmt::format("{}/{}", dir.string(), "mat.bin");

	if (!boost::filesystem::exists(mat_path)) { // write csr only if not exists already
		Eigen::write_binary_sparse<Eigen::SparseMatrix<T>>(mat_path, A);
		// return;
	}

	Eigen::write_binary_sparse<Eigen::SparseMatrix<T>>(mat_path, A);


	fmt::print("recip_edges: {}\n", recip_edges);
	fmt::print("recip_edges / m: {}\n", A.nonZeros());

	double reciprocity = double(recip_edges) / A.nonZeros();

	fmt::print("graph_name: {}\n", graph_name);
	fmt::print("sqlite3_db_path: {}\n", sqlite3_db_path);
	fmt::print("reciprocity: {}\n", reciprocity);
	single_val_set<double>(sqlite3_db_path, "n_vertices", "features", graph_name, n);
	single_val_set<double>(sqlite3_db_path, "n_edges", "features", graph_name, A.nonZeros());
	single_val_set<double>(sqlite3_db_path, "reciprocity", "features", graph_name, reciprocity);
// 	fmt::print("Solving directed..\n");
// 	double cyclic_eval = solve_directed(A, 1);
// 	fmt::print("Solving SVD..\n");
// 	double op_2_norm = compute_svd(A, 1);

// 	fmt::print("Symmetrizing..\n");
// 	// symmetrize the adjacencies
// 	A += Eigen::SparseMatrix<T>(A.transpose());
// #pragma omp parallel for schedule(dynamic)
// 	for (int k = 0; k < A.outerSize(); ++k)
// 		for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, k); it; ++it)
// 			if (it.valueRef() == 2) {
// 				it.valueRef() = 1;
// 			}

// 	fmt::print("Symmetrized. Solving Symmetric..\n");

// 	std::pair<double, double> p = solve_symmetric(A, 2);
// 	double spec_norm = p.first;
// 	double spec_sep = p.second;

// 	fmt::print("cyclic_eval: {}\n", cyclic_eval);
// 	fmt::print("op_2_norm: {}\n", op_2_norm);
// 	fmt::print("spec_norm: {}\n", spec_norm);
// 	fmt::print("spec_sep: {}\n", spec_sep);
	return;

}




//	single_val_set<double>(sqlite3_db_path, "op_2_norm", "features", graph_name, op_2_norm);
//	single_val_set<double>(sqlite3_db_path, "cyclic_eval", "features", graph_name, cyclic_eval);
//	single_val_set<double>(sqlite3_db_path, "spectral_norm", "features", graph_name, spec_norm);
//	single_val_set<double>(sqlite3_db_path, "spectral_separation", "features", graph_name, spec_sep);



#endif //GRAPH_PREPROCESS_STATS_H
