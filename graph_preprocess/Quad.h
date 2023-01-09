//
// Created by atrostan on 04/01/23.
//

#ifndef GRAPH_PREPROCESS_QUAD_H
#define GRAPH_PREPROCESS_QUAD_H
#include <cstdint>

class Quad {
public:
	uint32_t qx;
	uint32_t qy;
	uint32_t q_idx;
	uint32_t nnz;
	uint32_t *edges = nullptr;

	Quad() {
		nnz = 0;
	}

	Quad(Quad &&other) : qx(other.qx), qy(other.qy),
	                     q_idx(other.q_idx), nnz(other.nnz), edges(other.edges) {
		other.qx = -1;
		other.qy = -1;
		other.q_idx = -1;
		other.nnz = -1;
		other.edges = nullptr;
	}

	Quad &operator=(Quad &&other) {
		if (this != &other) {
			ReleaseResources();
			qx = other.qx;
			qy = other.qy;
			q_idx = other.q_idx;
			nnz = other.nnz;
			edges = other.edges;
			other.qx = -1;
			other.qy = -1;
			other.q_idx = -1;
			other.nnz = -1;
			other.edges = nullptr;
		}
		return *this;
	}

	void ReleaseResources() {
		if (edges != nullptr) {
			delete[] edges;
		}
	}

	~Quad() {
		ReleaseResources();
	}
};


#endif //GRAPH_PREPROCESS_QUAD_H
