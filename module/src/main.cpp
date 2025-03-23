#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include "smo.h"

namespace py = pybind11;


PYBIND11_MODULE(qpproblem_smo_solver, m)
{
    m.doc() = "C++ module solving QP problem appearing while training Support Vector Machine (binary classification) with regularization parametr C";

    py::class_<QPSolver>(m, "QPSolver")
        .def(py::init<
                 std::vector<std::vector<double>>,
                 std::vector<double>,
                 std::string,
                 double,
                 double,
                 int>(),
             py::arg("X"),
             py::arg("y"),
             py::arg("kernel"),
             py::arg("C"),
             py::arg("tol"),
             py::arg("max_iter"))
        .def("solve", &QPSolver::solve)
        .def("get_alpha", &QPSolver::get_alpha)
        .def("get_b", &QPSolver::get_b);
}