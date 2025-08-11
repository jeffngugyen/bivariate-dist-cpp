#ifndef DIST_H
#define DIST_H

#include <omp.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include <chrono>
#include <vector>
#include <cstdint>
#include <set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <bitset>
#include <cassert>
#include <compare>

std::tuple<int, int> compute_dist(std::vector<std::vector<int> > const& X, std::vector<std::vector<int> > const& Z, int l=4, int m=7);

#endif
