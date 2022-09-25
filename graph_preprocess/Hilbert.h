//
// Created by atrostan on 12/09/22.
//

#ifndef GRAPH_PREPROCESS_HILBERT_H
#define GRAPH_PREPROCESS_HILBERT_H

#include <cstdint>

void rot(int64_t n, int64_t* x, int64_t* y, int64_t rx, int64_t ry);

// convert (x,y) to d
int64_t xy2d(int64_t n, int64_t x, int64_t y);

// convert d to (x,y)
void d2xy(int64_t n, int64_t d, int64_t* x, int64_t* y);

int hyperceiling(int N);

#endif //GRAPH_PREPROCESS_HILBERT_H
