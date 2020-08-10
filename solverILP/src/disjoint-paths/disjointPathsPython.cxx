/*
 * disjointPathsPython.cxx
 *
 *  Created on: Jul 7, 2020
 *      Author: fuksova
 */

#include <pybind11/pybind11.h>
#include "disjoint-paths/disjointParams.hxx"
#include "disjoint-paths/disjointPathsMethods.hxx"
#include "disjoint-paths/ilp/solver-disjoint-ilp.hxx"

namespace py = pybind11;



PYBIND11_MODULE(disjointPathsPy, m) {
    m.doc() = "python binding for lifted disjoint paths";

     py::class_<disjointPaths::DisjointParams<>>(m, "DisjointParams")
        .def(py::init<const std::string &>());
        
     py::class_<disjointPaths::VertexGroups<>>(m, "TimeFramesToVertices")
        .def(py::init<>())
        .def("init_from_vector", &disjointPaths::VertexGroups<>::initFromVector, "Initializes vertices in time frames from a vector of size_t")
        .def("init_from_file", &disjointPaths::VertexGroups<>::initFromFile, "Initializes vertices in time frames from a file");

     py::class_<disjointPaths::CompleteStructure<>>(m, "GraphStructure")
        .def(py::init<disjointPaths::VertexGroups<> &>())
        .def("add_edges_from_array", &disjointPaths::CompleteStructure<>::addEdgesFromMatrix, "Initializes edges of the graph between two time frames from a matrix.")
        .def("add_edges_from_file", &disjointPaths::CompleteStructure<>::addEdgesFromFile, "Initializes all edges of the graph from a file.");

    // m.def("solve_ilp", py::overload_cast<disjointPaths::DisjointParams<>&, disjointPaths::CompleteStructure<>&>(&disjointPaths::solver_ilp_intervals<>), "Solve lifted disjoint paths");

     m.def("solve_ilp", &disjointPaths::solver_ilp<>, "Solve lifted disjoint paths");

     m.def("solve_standard_disjoit_paths", &disjointPaths::solver_flow_only<>, "Solve standard disjoint paths problem");




}


