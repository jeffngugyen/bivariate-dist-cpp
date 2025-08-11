#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "dist.h"

#include <tuple>
#include <vector>

namespace py = pybind11;

PYBIND11_MODULE(dist_py, m) {
    m.doc() = "my pybind11 plugin for computing hamming CSS matrix distances";

    m.def("compute_dist", &compute_dist, "Compute the hamming CSS matrix distance");
}
