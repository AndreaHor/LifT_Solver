#ifndef COMPLETESTRUCTURE_HXX
#define COMPLETESTRUCTURE_HXX

#include <stdexcept>
#include "disjoint-paths/disjointPathsMethods.hxx"
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <andres/graph/digraph.hxx>
#include <iterator>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <utility>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

namespace disjointPaths {



template<class T = size_t>
struct CompleteStructure {
public:
    template<class PAR>
    CompleteStructure(PAR & configParameters)
{
    //    std::cout<<"time file name check 1 "<<configParameters.getTimeFileName()<<std::endl;
        pVertexGroups=new VertexGroups<>(configParameters);
        VertexGroups<>& vg=*pVertexGroups;
                maxTime=vg.getMaxTime();

        configParameters.getControlOutput()<<"max time "<<maxTime<<std::endl;
                completeGraph=andres::graph::Digraph<>(vg.getMaxVertex()+1);

        configParameters.getControlOutput()<<"cg vertices "<<completeGraph.numberOfVertices()<<std::endl;
        addEdgesFromFile(configParameters.getGraphFileName(),configParameters);
        deleteVG=true;

}


    CompleteStructure(VertexGroups<>& vertexGroups):
        pVertexGroups(&vertexGroups)
{
        VertexGroups<>& vg=*pVertexGroups;
        maxTime=vg.getMaxTime();
        std::cout<<"max time "<<maxTime<<std::endl;
        completeGraph=andres::graph::Digraph<>(vg.getMaxVertex()+1);
        std::cout<<"cg vertices "<<completeGraph.numberOfVertices()<<std::endl;
        deleteVG=false;


}
    ~CompleteStructure(){
        if(deleteVG){
            delete pVertexGroups;
        }
    }


    template<class PAR>
    void addEdgesFromFile(const std::string& fileName,PAR& params);
    void addEdgesFromMatrix(size_t time1,size_t time2,const py::array_t<double> inputMatrix);
    template<class PAR>
    void addEdgesFromVectors(const py::array_t<size_t>edges,const py::array_t<double> costs,PAR& params);
    const VertexGroups<>& getVertexGroups(){
        return *pVertexGroups;
    }
    //void addEdges(size_t time1,size_t time2,std::stringstream& data);

        andres::graph::Digraph<> completeGraph;
        std::vector<double> completeScore;

        size_t maxTime;
    //DisjointParams<T>& params;

private:
    bool deleteVG;
    VertexGroups<>* pVertexGroups;

};

template<class T>
template<class PAR>
inline void CompleteStructure<T>::addEdgesFromVectors(const py::array_t<size_t> edges,const py::array_t<double> costs,PAR& params){
    char delim=',';
    VertexGroups<>& vg=*pVertexGroups;

    const auto edgeVector=edges.unchecked<2>();
    const std::size_t dim1=edgeVector.shape(0);
    const std::size_t dim2=edgeVector.shape(1);
    const auto costVector=costs.unchecked<1>();
    const size_t dimOfCosts=costVector.shape(0);


    if(dim2!=2){
        std::string message="Wrong dimension of edge array, second dimension 2 expected";
        throw std::invalid_argument(message);
    }
    if(dim1!=dimOfCosts){
        std::string message="Dimension of edge array and edge costs do not match.";
        throw std::invalid_argument(message);
    }

    params.getControlOutput()<<"Reading base edges from vector. "<<std::endl;
    params.writeControlOutput();
    for (size_t i=0;i<dim1;i++) {
        size_t v=edgeVector(i,0);
        size_t w=edgeVector(i,1);
        double edgeCost=costVector(i);

        size_t l0=vg.getGroupIndex(v);
        size_t l1=vg.getGroupIndex(w);
        //std::cout<<std::to_string(l0)<<", "<<std::to_string(l1)<<std::endl;
        if(v>vg.getMaxVertex()||w>vg.getMaxVertex()) continue;
        //std::cout<<"edge "<<v<<" "<<w<<std::endl;

        if(l1-l0<=params.getMaxTimeGapComplete()){
            completeGraph.insertEdge(v,w);
            completeScore.push_back(edgeCost);
        }
    }
}


template<class T>
inline void CompleteStructure<T>::addEdgesFromMatrix(size_t time1,size_t time2,const py::array_t<double> inputMatrix){

    VertexGroups<>& vg=*pVertexGroups;
    const auto matrix=inputMatrix.unchecked<2>();
    const std::size_t dim1=matrix.shape(0);
    const std::size_t dim2=matrix.shape(1);

    size_t transformIndex1=vg.getGroupVertices(time1)[0];
    size_t numberOfVertices1=vg.getGroupVertices(time1).size();
    size_t transformIndex2=vg.getGroupVertices(time2)[0];
    size_t numberOfVertices2=vg.getGroupVertices(time2).size();

    if(dim1!=numberOfVertices1||dim2!=numberOfVertices2){
        std::string message="Dimension mismatch, expected dimensions: "+std::to_string(numberOfVertices1)+", "+std::to_string(numberOfVertices2)+", got: "+std::to_string(dim1)+", "+std::to_string(dim2);
        throw std::invalid_argument(message);
    }
    for(std::size_t i=0; i<dim1; ++i) {
        size_t vertex1=i+transformIndex1;
        for(std::size_t j=0; j<dim2; ++j) {
            double score=matrix(i,j);
            if(!isinf(score)){
                size_t vertex2=j+transformIndex2;
                completeGraph.insertEdge(vertex1,vertex2);
                completeScore.push_back(score);
            }
        }
    }

}


//template<class T>
//inline void CompleteStructure<T>::addEdges(size_t time1,size_t time2,std::stringstream& data){
//	std::string line;
//	std::vector<std::string> strings;
//	char delim=',';

//	size_t lineCounter=0;
//	size_t transformIndex1=vg.getGroupVertices(time1)[0];
//	size_t numberOfVertices1=vg.getGroupVertices(time1).size();
//	size_t transformIndex2=vg.getGroupVertices(time2)[0];
//	size_t numberOfVertices2=vg.getGroupVertices(time2).size();

//	while (std::getline(data, line) && !line.empty()) {
//		if(lineCounter>=numberOfVertices1){
//			throw std::invalid_argument("Attempt to add edges to non-existing vertices in time frame "+time1);
//		}
//		size_t vertex1=lineCounter+transformIndex1;
//		strings = split(line, delim);
//		if(strings.size()!=numberOfVertices2){
//			throw std::invalid_argument("Number of matrix columns does not match the number of vertices in time frame "+time2);
//		}

//		for (int i = 0; i < numberOfVertices2; ++i) {
//			double score=std::stoul(strings[i]);
//			if(score<std::numeric_limits<double>::infinity()){
//				size_t vertex2=i+transformIndex2;
//				completeGraph.insertEdge(vertex1,vertex2);
//				completeScore.push_back(score);
//			}
//		}
//    	lineCounter++;
//	}
//	if(lineCounter<numberOfVertices1){
//		throw std::invalid_argument("Number of lines in matrix is lower than number of vertices in time "+time1);
//	}
//}


template<class T>
template<class PAR>
inline void CompleteStructure<T>::addEdgesFromFile(const std::string& fileName,PAR& params){
        char delim=',';
    VertexGroups<>& vg=*pVertexGroups;

        std::string line;
        std::ifstream data;
        try{
        data.open(fileName);
                if(!data){
            throw std::system_error(errno, std::system_category(), "failed to open graph file "+fileName);
                }

                std::getline(data, line);
                double objValue=0;


        params.getControlOutput()<<"Read big graph" << std::endl;
                std::vector<std::string> strings;

        params.getControlOutput()<<"Reading vertices from file. "<<std::endl;
        params.writeControlOutput();
                //Vertices that are not found have score=0. Appearance and disappearance cost are read here.
                while (std::getline(data, line) && !line.empty()) {

                }

        params.getControlOutput()<<"Reading base edges from file. "<<std::endl;
        params.writeControlOutput();
                size_t maxLabel=0;
                while (std::getline(data, line) && !line.empty()) {

            //std::cout<<line<<std::endl;
                        strings = split(line, delim);

                        unsigned int v = std::stoul(strings[0]);
            //std::cout<<v<<std::endl;
                        unsigned int w = std::stoul(strings[1]);
            //std::cout<<w<<std::endl;
                        size_t l0=vg.getGroupIndex(v);
                        size_t l1=vg.getGroupIndex(w);
            //std::cout<<std::to_string(l0)<<", "<<std::to_string(l1)<<std::endl;
                        if(v>vg.getMaxVertex()||w>vg.getMaxVertex()) continue;
                        //std::cout<<"edge "<<v<<" "<<w<<std::endl;

                        if(l1-l0<=params.getMaxTimeGapComplete()){
                                double score = std::stod(strings[2]);
                                completeGraph.insertEdge(v,w);
                                completeScore.push_back(score);
                        }

                }


                //objValue+=maxLabel*parameters.getInputCost()+maxLabel*parameters.getOutputCost();

                data.close();
        }
        catch (std::system_error& er) {
                std::clog << er.what() << " (" << er.code() << ")" << std::endl;

        }
}

}//End of namespace
#endif // COMPLETESTRUCTURE_HXX
