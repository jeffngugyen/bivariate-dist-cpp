#include "dist.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(dist_py, module) {
    module.doc() = "Exact small CSS-code distances using binary linear algebra";
    module.def("compute_dist", &compute_dist,
               py::arg("X"), py::arg("Z"), py::arg("l") = 4, py::arg("m") = 7,
               py::call_guard<py::gil_scoped_release>(),
               "Return (d_Z, d_X) for X/Z stabilizer matrices with 2*l*m columns. "
               "Preserves the original tuple order. Invalid inputs raise ValueError.");
}
