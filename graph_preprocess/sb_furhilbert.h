//
// Created by atrostan on 14/09/22.
//

#include <iostream>
//#include "furhilbert.h"
#include "typedefs.h"
#include <chrono>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <unistd.h>
#include <fstream>
#include <omp.h>
#include <vector>
#include "QuadrantCalculator.h"
#include "util.h"

#ifndef GRAPH_PREPROCESS_SB_FURHILBERT_H
#define GRAPH_PREPROCESS_SB_FURHILBERT_H

void compute_quadrants(uint32_t n);

void compute_section(int depth, int rot, uint32_t start_x, uint32_t end_x,
                     uint32_t start_y, uint32_t end_y, int critical_depth,
                     std::vector<Quadrant> &qs, int qidx);

#endif //GRAPH_PREPROCESS_SB_FURHILBERT_H
