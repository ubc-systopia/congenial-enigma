//
// Created by atrostan on 28/12/22.
//
#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/numeric>

#include <getopt.h>
#include <cstdint>
#include <vector>
#include <string>
#include <omp.h>
#include <likwid-marker.h>
#include <likwid.h>
#include "fmt/core.h"
#include "fmt/ranges.h"
#include "io.h"
#include "omp.h"
#include "util.h"
#include <bitset>
#include <limits.h>
#include "ips4o/ips4o.hpp"
#include "Quad.h"
#include "CSR.h"

#ifndef GRAPH_PREPROCESS_HILBURN_PREPROC_H
#define GRAPH_PREPROCESS_HILBURN_PREPROC_H


//std::pair<uint64_t, uint64_t> bsearch(uint32_t *arr, uint64_t x, uint64_t l, uint64_t h) {
//	uint64_t mid;
//	while (h - l > 1) {
//		mid = (h + l) / 2;
//		if (arr[mid] < x) l = mid + 1;
//		else h = mid;
//	}
//	return {l, h};
////	if (arr[l] == x) return l;
////	else if (arr[h] == x) return h;
////	else return -1;
//}



#endif //GRAPH_PREPROCESS_HILBURN_PREPROC_H
