#ifndef TWO_GRAPHS_INPUT_STRUCTURE_HXX_
#define TWO_GRAPHS_INPUT_STRUCTURE_HXX_

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include "disjoint-paths/vertexGroups.hxx"

namespace py = pybind11;
namespace disjointPaths {

class TwoGraphsInputStructure
{
public:
    TwoGraphsInputStructure(const py::array_t<size_t>& baseEdges,const py::array_t<size_t>& liftedEdges,const  py::array_t<double>& baseCosts,const  py::array_t<double>& liftedCosts,const VertexGroups<>* pvg) {
        pBaseEdges=&baseEdges;
        pLiftedEdges=&liftedEdges;
        pBaseCosts=&baseCosts;
        pLiftedCosts=&liftedCosts;
        myPvg=pvg;
    }
    const py::array_t<size_t>* pBaseEdges;
    const py::array_t<size_t>* pLiftedEdges;
    const py::array_t<double>* pBaseCosts;
    const py::array_t<double>* pLiftedCosts;
    const VertexGroups<>* myPvg;


};

}

#endif // TWOGRAPHSINPUTSTRUCTURE_HXX_
