/*
 * disjointStructure.hxx
 *
 *  Created on: Jul 3, 2020
 *      Author: fuksova
 */

#ifndef INCLUDE_DISJOINT_PATHS_DISJOINTSTRUCTURE_HXX_
#define INCLUDE_DISJOINT_PATHS_DISJOINTSTRUCTURE_HXX_


#include <stdexcept>
#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <andres/graph/digraph.hxx>
#include <stack>
#include <unordered_set>
#include <iterator>
#include <unordered_map>
#include <string>
#include <bitset>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "levinkov/timer.hxx"
#include <deque>
#include <queue>
#include <set>
#include <map>
#include <list>
#include "disjoint-paths/disjointParams.hxx"
#include "disjoint-paths/disjointPathsMethods.hxx"
#include "disjoint-paths/completeStructure.hxx"
#include <utility>

namespace disjointPaths {

template<class T = size_t>
struct DisjointStructure {
public:
	typedef T size_type;


	DisjointStructure(disjointPaths::DisjointParams<>& configParameters,char delim=',',CompleteStructure<>* cs=0,size_t minTime=0,size_t maxTime=0);
//	bool isReachable(size_t i,size_t j) const{
//		if(i==t_) return false;  //Assume no path from the terminal node
//		if(desc.size()==0) return true;
//		return desc[i][j];


//	}

    bool isReachable(size_t i,size_t j) const{
		if(i==t_) return false;  //Assume no path from the terminal node
		if(reachable.size()==0) return true;
		return reachable[i].count(j)>0;
	}

	std::unordered_set<size_t>& reachableFromVertex(size_t v){
		return reachable[v];
	}

//	std::vector<std::vector<bool>> getReachable(){
//		return desc;
//	}

	andres::graph::Digraph<>& getGraph()  {
		return graph_;
	}
	andres::graph::Digraph<>& getGraphLifted()  {
		return graphLifted_;
	}


    void setGraph(const andres::graph::Digraph<>& newGraph)  {
        graph_=newGraph;
    }

    void setGraphLifted(const andres::graph::Digraph<>& newGraphLifted)  {
        graphLifted_=newGraphLifted;
    }

//	std::vector<std::vector<bool>>* getPReachable(){
//		return &desc;
//	}

    std::vector<std::unordered_set<size_t>>* getPReachable(){
		return &reachable;
	}

	andres::graph::Digraph<>* getPGraph()  {
		return &graph_;
	}
	andres::graph::Digraph<>* getPGraphLifted()  {
		return &graphLifted_;
	}

//	VertexGroups<>& getVertexGroups(){
//		return vertexGroups;
//	}

	size_t getSourceNode() const {
		return s_;
	}

	size_t getTerminalNode() const {
		return t_;
	}

	double getEdgeScore(size_t e) const {
		return edgeScore[e];
	}

	double getEdgeScore(size_t v0,size_t v1) const {
		auto findEdge=graph_.findEdge(v0,v1);
		assert(findEdge.first);
		return edgeScore[findEdge.second];
	}

	const std::vector<double>& getEdgesScore() {
		return edgeScore;
	}

    void setEdgesScore(const std::vector<double>& beScore){
        edgeScore=beScore;
    }


	const std::vector<double>& getLiftedEdgesScore() {
			return liftedEdgeScore;
	}

    void setLiftedEdgesScore(const std::vector<double>& leScore){
        liftedEdgeScore=leScore;
    }


	const std::vector<double>& getVerticesScore() {
		return vertexScore;
	}

	double getLiftedEdgeScore(size_t e) const {
		return liftedEdgeScore[e];
	}



	double getLiftedEdgeScore(size_t v0,size_t v1) const {
		auto findEdge=graphLifted_.findEdge(v0,v1);
		assert(findEdge.first);
		return liftedEdgeScore[findEdge.second];
	}

	double getVertexScore(size_t v) const {
		return vertexScore[v];
	}

	VertexGroups<size_t>* getPVertexGroups(){
		return &vertexGroups;
	}


	disjointPaths::DisjointParams<>& parameters;

	bool isTimeFramesEnabled(){
		return useTimeFrames;
	}




    const VertexGroups<size_t>& getVertexGroups() const{
        return vertexGroups;
    }

    const size_t& getMinV() const{
        return minV;
    }


    const size_t& getMaxV() const{
        return maxV;
    }

private:

    size_t s_;
    size_t t_;

    VertexGroups<size_t> vertexGroups;
    size_t minV=0;
    size_t maxV=0;

	std::vector<double> vertexScore;
	std::vector<double> edgeScore;
	std::vector<double> liftedEdgeScore;
    //std::vector<std::vector<bool>> desc;  //TODO remove
	std::vector<std::unordered_set<size_t>> reachable;



	andres::graph::Digraph<> graph_;
	andres::graph::Digraph<> graphLifted_;


	void readGraphWithTime(size_t minTime,size_t maxTime,CompleteStructure<>* cs);

	bool useTimeFrames;


};


template<class T>
inline DisjointStructure<T>::DisjointStructure(disjointPaths::DisjointParams<>& configParameters,char delim,CompleteStructure<>* cs,size_t minTime,size_t maxTime):
parameters(configParameters)
{

    parameters.getControlOutput()<<"interval "<<minTime<<","<<maxTime<<std::endl;
    parameters.writeControlOutput();

	useTimeFrames=parameters.isRestrictFrames()||parameters.isSparsify();
	size_t maxVertex;


    readGraphWithTime(minTime,maxTime,cs);



	if(graph_.numberOfVertices()>2){

		if(parameters.getMaxTimeLifted()>0&&parameters.isAutomaticLifted()){

            parameters.getControlOutput()<<"Adding automatic lifted edges"<<std::endl;
            parameters.writeControlOutput();
			for (int i = 0; i < graph_.numberOfEdges(); ++i) {
				size_t v0=graph_.vertexOfEdge(i,0);
				size_t v1=graph_.vertexOfEdge(i,1);
				if(v0!=s_&&v1!=t_){
					//	if(secOrderDesc[v0][v1]){
					graphLifted_.insertEdge(v0,v1);
					liftedEdgeScore.push_back(edgeScore[i]);
					//				}
					//				else{
					//					edgeScore[i]*=2;
					//				}
				}
			}

            parameters.getControlOutput()<<"done"<<std::endl;
            parameters.writeControlOutput();

		}


        parameters.getControlOutput()<<"number of vertices "<<graph_.numberOfVertices()<<std::endl;
        parameters.writeControlOutput();
		if(parameters.isSparsify()){
            createKnnBaseGraph(*this,parameters);


            if(parameters.getMaxTimeLifted()>0){
                if(parameters.getSmallIntervals()==0){
                    //desc=initReachable(graph_,parameters,&vertexGroups);
                    reachable=initReachableSet(graph_,parameters,&vertexGroups);
                }
                else{
                    reachable=initReachableSet(graph_,parameters);
                }
                 keepFractionOfLifted(*this,parameters);

            }
            else{
                parameters.getControlOutput()<<"FW skipped"<<std::endl;
                parameters.writeControlOutput();
            }

    }
		else{
			//desc=initReachable(graph_,parameters);
			reachable=initReachableSet(graph_,parameters);
            //TODO remove redundant edges from lifted graph
		}
	}
}



//includes minTime, excludes maxTime
template<class T>
inline void DisjointStructure<T>::readGraphWithTime(size_t minTime,size_t maxTime,CompleteStructure<>* cs){

	andres::graph::Digraph<>& completeGraph=cs->completeGraph;
	std::vector<double>& completeScore=cs->completeScore;
    const VertexGroups<>& vg=cs->getVertexGroups();

	std::unordered_map<size_t,std::vector<size_t>> groups;

	size_t mt=minTime;
	while(vg.getGroupVertices(mt).size()==0&&mt<maxTime){
		mt++;
	}

	size_t minVertexTime=mt;

	if(minVertexTime==maxTime){
		//TODO all empty, do something meaningful
		s_ =0;
		t_ = 1;
		graphLifted_ = andres::graph::Digraph<>(2);
		graph_ = andres::graph::Digraph<>(2);
        parameters.getControlOutput()<<"empty interval in ds constructor"<<std::endl;
        parameters.writeControlOutput();

	}
	else{
        parameters.getControlOutput()<<"valid interval in ds constructor"<<std::endl;
        parameters.writeControlOutput();
		size_t minVertex=vg.getGroupVertices(minVertexTime)[0];

		mt=maxTime-1;
		while(vg.getGroupVertices(mt).size()==0&&mt>minVertexTime){
			mt--;
		}
		size_t maxVertexTime=mt;
		//		if(maxVertexTime==minVertexTime){
		//			//TODO graphs with one layer and no edges, do something meaningful
		//		}
		//		else{

		size_t maxVertex=*(vg.getGroupVertices(maxVertexTime).rbegin());

		size_t numberOfVertices=maxVertex-minVertex+3;
		s_ = numberOfVertices - 2;
		t_ = numberOfVertices - 1;

		std::vector<size_t> vToGroup(numberOfVertices);
		vToGroup[s_]=0;
		vToGroup[t_]=maxTime-minTime+1;
        std::vector<size_t> startGroup(1);
        startGroup[0]=s_;
        std::vector<size_t> terminalGroup(1);
        terminalGroup[0]=t_;
        groups[0]=startGroup;
        groups[maxTime-minTime+1]=terminalGroup;


		for (int gi = minTime; gi < maxTime; ++gi) {
			//groups[gi-minTime+1]=std::vector<size_t>();
			for(size_t v:vg.getGroupVertices(gi)){
				size_t vertex=v-minVertex;
				groups[gi-minTime+1].push_back(vertex);
				vToGroup[vertex]=gi-minTime+1;
			}
		}

		vertexGroups=VertexGroups<>(groups,vToGroup);
		vertexScore = std::vector<double>(numberOfVertices, 0);
		graphLifted_ = andres::graph::Digraph<>(numberOfVertices);
		graph_ = andres::graph::Digraph<>(numberOfVertices);

		bool useZeroInOut=false;
		for (int v = 0; v < numberOfVertices-2; ++v) {
			graph_.insertEdge(s_,v);
			edgeScore.push_back(parameters.getInputCost());
			graph_.insertEdge(v,t_);
			edgeScore.push_back(parameters.getOutputCost());
		}


		for (int v = minVertex; v < maxVertex; ++v) {
			for (int i = 0; i < completeGraph.numberOfEdgesFromVertex(v); ++i) {
				size_t w=completeGraph.vertexFromVertex(v,i);
				if(w>maxVertex) continue;
				size_t e=completeGraph.edgeFromVertex(v,i);
				graph_.insertEdge(v-minVertex, w-minVertex);
				edgeScore.push_back(completeScore[e]);
			}
		}
		minV=minVertex;
		maxV=maxVertex;
        parameters.getControlOutput()<<"DS constructor max Vertex "<<maxV<<std::endl;
        parameters.writeControlOutput();
		//}
	}

}



}





#endif /* INCLUDE_DISJOINT_PATHS_DISJOINTSTRUCTURE_HXX_ */
