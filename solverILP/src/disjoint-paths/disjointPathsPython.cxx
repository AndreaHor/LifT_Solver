/*
 * disjointPathsPython.cxx
 *
 *  Created on: Jul 7, 2020
 *      Author: fuksova
 */

#include <pybind11/pybind11.h>
#include "disjoint-paths/disjointParams.hxx"

namespace py = pybind11;



PYBIND11_MODULE(disjointPathsPy, m) {
    m.doc() = "python binding for lifted disjoint paths";

     py::class_<disjointPaths::DisjointParams<>>(m, "DisjointParams")
        .def(py::init<const std::string &>());
        

}


