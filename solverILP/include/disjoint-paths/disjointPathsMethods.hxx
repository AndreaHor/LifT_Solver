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




}


#endif /* INCLUDE_DISJOINT_PATHS_DISJOINTPATHSMETHODS_HXX_ */
