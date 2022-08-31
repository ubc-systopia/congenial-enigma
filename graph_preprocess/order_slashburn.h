//
// Created by atrostan on 30/08/22.
//

#include <vector>
#include "typedefs.h"
#include "igraph.h"
#include "util.h"

#ifndef GRAPH_PREPROCESS_ORDER_SLASHBURN_H
#define GRAPH_PREPROCESS_ORDER_SLASHBURN_H


std::tuple<igraph_t *, ul, ul> order_igraph_slashburn(igraph_t &g, const int k, std::vector<ul> &rank, ul &hub_idx, ul &spokes_end_idx);

#endif //GRAPH_PREPROCESS_ORDER_SLASHBURN_H


