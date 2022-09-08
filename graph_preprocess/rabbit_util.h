//
// Created by atrostan on 07/09/22.
//
#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <boost/range/algorithm/copy.hpp>
#include <iostream>
#include <deque>
#include "io.h"

#ifndef GRAPH_PREPROCESS_RABBIT_UTIL_H
#define GRAPH_PREPROCESS_RABBIT_UTIL_H


//typedef std::tuple<ul, ul, float> edge;
//typedef std::tuple<ul, ul> edge;
#define RO_DIE rabbit_order::aux::die_t(__FILE__, __LINE__, __func__)
std::vector<std::pair<ul, ul>> par_read(const std::string &path);

#endif //GRAPH_PREPROCESS_RABBIT_UTIL_H
