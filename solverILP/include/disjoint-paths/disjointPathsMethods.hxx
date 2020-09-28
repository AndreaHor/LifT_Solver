/*
 * disjointPathsMethods.hxx
 *
 *  Created on: Jul 3, 2020
 *      Author: fuksova
 */

#ifndef INCLUDE_DISJOINT_PATHS_DISJOINTPATHSMETHODS_HXX_
#define INCLUDE_DISJOINT_PATHS_DISJOINTPATHSMETHODS_HXX_



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
//#include "disjoint-paths/disjointParams.hxx"
#include <utility>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include "disjoint-paths/vertexGroups.hxx"

namespace py = pybind11;
namespace disjointPaths {



template<class T,class PAR>
    std::vector<std::unordered_set<size_t>> initReachableSet(T & graph,PAR& parameters,VertexGroups<size_t>* vg=0){

	levinkov::Timer tfw;
			tfw.start();

        parameters.getControlOutput()<<"Run Floyd Warshall"<<std::endl;
		const size_t n=graph.numberOfVertices();

		// todo: use some matrix structure
		std::vector<std::unordered_set<size_t>> desc(n);
		std::vector<std::vector<std::bitset<10000>>> descBit(n);
		size_t columns=n/10000;
		size_t mod=n%10000;
		if(mod>0) columns++;


		for (int i = 0; i < n; ++i) {
			descBit[i]=std::vector<std::bitset<10000>>(columns);
		}

		for (size_t v = 0; v < n; ++v) {

			descBit[v][v/10000][v%10000]=1; //make this reflexive
			size_t edges=graph.numberOfEdgesFromVertex(v);
			for (int j = 0; j < edges; ++j) {
				size_t w=graph.vertexFromVertex(v,j);
				descBit[v][w/10000][w%10000]=1;


			}
		}



		if(vg==0){
			for (int k1 = 0; k1 <columns; ++k1) {
				for (int k2 = 0; k2 < 10000; ++k2) {
					if(k1*10000+k2<n){
						for (int i = 0; i < n; ++i) {
							if(descBit[i][k1][k2]){
								for (int j = 0; j < columns; ++j) {
									descBit[i][j]|=descBit[k1*10000+k2][j];
								}

							}
						}
					}
					else{
						break;
					}
				}

			}
		}
		else{
			for (int k1 = 0; k1 <columns; ++k1) {
				for (int k2 = 0; k2 < 10000; ++k2) {
					if(k1*10000+k2<n){
						size_t maxTime=vg->getGroupIndex(k1*10000+k2);
						for (int t = 0; t < maxTime; ++t) {//TODO use time gap
							const std::vector<size_t>& vertices=vg->getGroupVertices(t);
							for (size_t i:vertices) {
								if(descBit[i][k1][k2]){
									for (int j = 0; j < columns; ++j) {
										descBit[i][j]|=descBit[k1*10000+k2][j];
									}

								}
							}
						}
					}
					else{
						break;
					}
				}

			}

		}

		for (int i = 0; i < n; ++i) {
			for (int k1 = 0; k1 <columns; ++k1) {
				for (int k2 = 0; k2 < 10000; ++k2) {
					if(k1*10000+k2<n&&descBit[i][k1][k2]){
						desc[i].insert(k1*10000+k2);
					}
				}
			}

		}

		tfw.stop();


        parameters.getControlOutput()<<"fw finished in time "<<tfw.get_elapsed_seconds() << std::endl;
		return desc;

}



template<class T,class PAR>
    std::vector<std::vector<bool>> initReachable(T & graph,PAR& parameters,VertexGroups<size_t>* vg=0){

	levinkov::Timer tfw;
			tfw.start();

        parameters.getControlOutput()<<"Run Floyd Warshall"<<std::endl;
		const size_t n=graph.numberOfVertices();

		// todo: use some matrix structure
		std::vector<std::vector<bool>> desc(n);
		std::vector<std::vector<std::bitset<10000>>> descBit(n);
		size_t columns=n/10000;
		size_t mod=n%10000;
		if(mod>0) columns++;


		for (int i = 0; i < n; ++i) {
			desc[i]=std::vector<bool>(n,0);
			descBit[i]=std::vector<std::bitset<10000>>(columns);
		}

		for (size_t v = 0; v < n; ++v) {

			descBit[v][v/10000][v%10000]=1; //make this reflexive
			size_t edges=graph.numberOfEdgesFromVertex(v);
			for (int j = 0; j < edges; ++j) {
				size_t w=graph.vertexFromVertex(v,j);
				descBit[v][w/10000][w%10000]=1;


			}
		}



		if(vg==0){
			for (int k1 = 0; k1 <columns; ++k1) {
				for (int k2 = 0; k2 < 10000; ++k2) {
					if(k1*10000+k2<n){
						for (int i = 0; i < n; ++i) {
							if(descBit[i][k1][k2]){
								for (int j = 0; j < columns; ++j) {
									descBit[i][j]|=descBit[k1*10000+k2][j];
								}

							}
						}
					}
					else{
						break;
					}
				}

			}
		}
		else{
			for (int k1 = 0; k1 <columns; ++k1) {
				for (int k2 = 0; k2 < 10000; ++k2) {
					if(k1*10000+k2<n){
						size_t maxTime=vg->getGroupIndex(k1*10000+k2);
						for (int t = 0; t < maxTime; ++t) {//TODO use time gap
							const std::vector<size_t>& vertices=vg->getGroupVertices(t);
							for (size_t i:vertices) {
								if(descBit[i][k1][k2]){
									for (int j = 0; j < columns; ++j) {
										descBit[i][j]|=descBit[k1*10000+k2][j];
									}

								}
							}
						}
					}
					else{
						break;
					}
				}

			}

		}

		for (int i = 0; i < n; ++i) {
			for (int k1 = 0; k1 <columns; ++k1) {
				for (int k2 = 0; k2 < 10000; ++k2) {

					if(k1*10000+k2<n){
						desc[i][k1*10000+k2]=descBit[i][k1][k2];
					}
				}
			}

		}

		tfw.stop();

        parameters.getControlOutput()<<"fw finished in time "<<tfw.get_elapsed_seconds() << std::endl;
		return desc;

}



//    template<class INSTANCE, class PAR>
//    inline void keepFractionOfLifted(INSTANCE& instance, const PAR& parameters){


//        parameters.getControlOutput()<<"Sparsify lifted graph"<<std::endl;
//        parameters.writeControlOutput();
//        //TODO run automaticLifted to find candidates first

//        double negMaxValue=0;
//        double posMinValue=0;
//        bool useAdaptive=false;
//        //TODO adaptive lifted threshold in config file
//        std::vector<double> newLiftedCosts;

//        const andres::graph::Digraph<>& graph_=instance.getGraph();
//        const andres::graph::Digraph<>& graphLifted_=instance.getGraphLifted();
//        const std::vector<double>& liftedCosts=instance.getLiftedEdgesScore;
//        const std::vector<std::unordered_set<size_t>>* pReachable =instance.getPReachableNew();
//        const std::vector<std::unordered_set<size_t>>& reachable=*pReachable;
//        const size_t t_=instance.getTerminalNode();

//        std::vector<double> newBaseEdgeScore=instance.getEdgesScore();


//        andres::graph::Digraph<> tempGraphLifted=(graph_.numberOfVertices());


//        if(useAdaptive){

//            double fPositive=parameters.getPositiveThresholdLifted();
//            double fNegative=parameters.getNegativeThresholdLifted();

//            assert(fNegative>=0&&fNegative<=1);
//            assert(fPositive>=0&&fPositive<=1);

//            std::priority_queue<double> allPositive;
//            std::priority_queue<double> allNegative;
//            for (int i = 0; i < graphLifted_.numberOfEdges(); ++i) {
//                double cost=liftedCosts.at(i);
//                if(cost<0){
//                    allNegative.push(-cost);
//                }
//                else{
//                    allPositive.push(cost);
//                }
//            }

//            double negativeToKeep=fNegative*allNegative.size();
//            double positiveToKeep=fPositive*allPositive.size();
//            int negNumber=int(round(negativeToKeep));
//            int posNumber=int(round(positiveToKeep));

//            for (int i = 0; i < negNumber; ++i) {
//                allNegative.pop();
//            }
//            for (int i = 0; i < posNumber; ++i) {
//                allPositive.pop();
//            }



//            if(!allNegative.empty()){
//                negMaxValue=-allNegative.top();
//            }

//            parameters.getControlOutput()<<"max negative lifted "<<negMaxValue<<std::endl;

//            if(!allPositive.empty()){
//                posMinValue=allPositive.top();
//            }

//            parameters.getControlOutput()<<"min positive lifted "<<posMinValue<<std::endl;
//            parameters.writeControlOutput();

//        }
//        else{
//            negMaxValue=parameters.getNegativeThresholdLifted();
//            posMinValue=parameters.getPositiveThresholdLifted();
//        }


//        std::unordered_map<size_t,std::set<size_t>> liftedEdges;
//        for (int v = 0; v < graphLifted_.numberOfVertices()-2; ++v) {
//            std::unordered_set<size_t> alternativePath;
//            for (int i = 0; i < graph_.numberOfEdgesFromVertex(v); ++i) {
//                size_t w=graph_.vertexFromVertex(v,i);
//                for(size_t u:reachable[w]){
//                    if(u!=w) alternativePath.insert(u);
//                }
//            }
//            for (int i = 0; i < graphLifted_.numberOfEdgesFromVertex(v); ++i) {
//                size_t w=graphLifted_.vertexFromVertex(v,i);
//                if(w!=t_){
//                    if(alternativePath.count(w)>0) liftedEdges[v].insert(w);

//                }

//            }
//        }


//        parameters.getControlOutput()<<"done"<<std::endl;
//        parameters.writeControlOutput();


//        for (int i = 0; i < graphLifted_.numberOfEdges(); ++i) {
//            size_t v0=graphLifted_.vertexOfEdge(i,0);
//            size_t v1=graphLifted_.vertexOfEdge(i,1);
//            int l0=vertexGroups.getGroupIndex(v0);
//            int l1=vertexGroups.getGroupIndex(v1);
//            double cost=liftedCosts.at(i);
//            bool goodCost=(cost<=negMaxValue)||(cost>=posMinValue);
//            //if(isReachable(v0,v1)){
//            if(instance.isReachableNew(v0,v1)){

//                int timeGapDiff=l1-l0-parameters.getDenseTimeLifted();
//                bool timeConstraint=l1-l0<=parameters.getDenseTimeLifted()||((l1-l0)<=parameters.getMaxTimeLifted()&&(timeGapDiff%parameters.getLongerIntervalLifted())==0);
//                //int modulo=timeGapDiff%parameters.getLongerIntervalLifted();
//                //if(l1-l0<=parameters.getDenseTimeLifted()||((l1-l0)<=parameters.getMaxTimeLifted()&&(timeGapDiff%parameters.getLongerIntervalLifted())==0)){
//                if(timeConstraint&&goodCost){
//                    if(liftedEdges[v0].count(v1)>0){
//                        tempGraphLifted.insertEdge(v0,v1);
//                        newLiftedCosts.push_back(cost);
//                    }
//                    else{
//                        auto edgeTest=graph_.findEdge(v0,v1);
//                        if(edgeTest.first){
//                            newBaseEdgeScore[edgeTest.second]+=cost;  //Compensate that the lifted edge has been removed
//                        }
//                        //						else{
//                        //							throw std::runtime_error("Error in lifted graph sparsification.");
//                        //						}
//                    }

//                }
//            }

//        }




////TODO set these via setters

//        edgeScore=newBaseEdgeScore;

//        liftedEdgeScore=newLiftedCosts;

//        graphLifted_=tempGraphLifted;

//        parameters.getControlOutput()<<"Left "<<newLiftedCosts.size()<<" lifted edges."<<std::endl;
//        parameters.writeControlOutput();

//        if(graphLifted_.numberOfEdges()!=newLiftedCosts.size()){

//            parameters.getControlOutput()<<"lifted edge number mismatch, lifted graph: "<<graphLifted_.numberOfEdges()<<", cost vector "<<newLiftedCosts.size()<<std::endl;
//        }
//        else{

//            parameters.getControlOutput()<<"lifted edge number and lifted graph size match "<<std::endl;

//        }
//        parameters.writeControlOutput();

//    }




}


#endif /* INCLUDE_DISJOINT_PATHS_DISJOINTPATHSMETHODS_HXX_ */
