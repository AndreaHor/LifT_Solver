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
#include "disjoint-paths/disjointParams.hxx"
#include <utility>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
namespace disjointPaths {



template<class T=size_t>
std::vector<std::vector<T>> readLines(std::string inputFileName, char delim) {
	std::ifstream inputFile;
	std::vector<std::vector<T>> outputs;
	try{

		inputFile.open(inputFileName);
		if (!inputFile){
			throw std::system_error(errno, std::system_category(), "failed to open "+inputFileName);
		}

		std::string line;
		size_t lineCounter=0;
		std::getline(inputFile, line);
		std::vector<std::string> strings;

		while (std::getline(inputFile, line) && !line.empty()) {
			strings=split(line,delim);
			size_t length=strings.size();
			std::vector<T> parsedLine(length);
			for (int i = 0; i < length; ++i) {
				parsedLine[i]=std::stoul(strings[i]);
			}
			outputs.push_back(parsedLine);

		}

		inputFile.close();


	}

	catch (std::system_error& er) {
		std::clog << er.what() << " (" << er.code() << ")" << std::endl;

	}
	return outputs;


}


template<class T=char>
 std::vector<std::string> split(
		std::string inputString, T delim) {
	size_t occurence = 0;
	size_t newOccurence = 0;
	std::vector<std::string> strings;
	while (newOccurence < inputString.size()) {
		newOccurence = std::min(inputString.find_first_of(delim, occurence),
				inputString.size());

		std::string newString(inputString, occurence, newOccurence - occurence);
		strings.push_back(newString);
		newOccurence = newOccurence + 1;
		occurence = newOccurence;
	}

	return strings;
}



template<class T = size_t>
struct VertexGroups {
public:
	VertexGroups(std::unordered_map<size_t,std::vector<size_t>> groups_,std::vector<size_t> vToGroup_):
		groups(groups_),vToGroup(vToGroup_)
	{
		maxVertex=vToGroup_.size()-3;
		maxTime=*(vToGroup.rbegin());
	}
	VertexGroups(){
		maxVertex=0;
		maxTime=0;
	}



	VertexGroups(disjointPaths::DisjointParams<>& parameters,char delim){
		//std::vector<std::vector<size_t>> groups;

		std::string line;
		std::vector<std::string> strings;
		std::ifstream timeData;
		try{
			timeData.open(parameters.getTimeFileName());
			if(!timeData){
				throw std::system_error(errno, std::system_category(), "failed to open "+parameters.getTimeFileName());
			}

			initFromStream(timeData,parameters.getMaxTimeFrame());

		}
		catch (std::system_error& er) {
			std::clog << er.what() << " (" << er.code() << ")" << std::endl;

		}

	}


	template<class STR>
	void initFromStream(STR& timeData,size_t maxTimeToRead);

    void initFromVector(const std::vector<size_t>& verticesInFrames);




	size_t getGroupIndex(size_t v){
		return vToGroup[v];
	}

	const std::vector<size_t>& getGroupVertices(size_t index){
		return groups[index];
	}

	//ID of maximal valid vertex (i.e. without s and t)
	size_t getMaxVertex() const {
		return maxVertex;
	}

	//Time of the last video frame
	size_t getMaxTime() const {
			return maxTime;
    }




private:
	//std::vector<std::vector<size_t>> groups;
	std::unordered_map<size_t,std::vector<size_t>> groups;
	std::vector<size_t> vToGroup;
	size_t maxVertex;
	size_t maxTime;


};


template<class T>
inline void VertexGroups<T>::initFromVector(const std::vector<size_t>& verticesInFrames){
    maxTime=verticesInFrames.size();
    size_t inFrameCounter=0;
    size_t vertexCounter=0;
    size_t frameCounter=1;
    std::vector<size_t> verticesInGroup;
    vToGroup=std::vector<size_t>();
    while(frameCounter<=maxTime){
        while(inFrameCounter==verticesInFrames[frameCounter-1]){
            groups[frameCounter]=verticesInGroup;
            inFrameCounter=0;
            frameCounter++;
            verticesInGroup=std::vector<size_t>();
        }
        if(frameCounter<=maxTime){
            verticesInGroup.push_back(vertexCounter);
            vToGroup.push_back(frameCounter);
            inFrameCounter++;
            vertexCounter++;
        }
    }

    maxVertex=vertexCounter;
    size_t s=maxVertex+1;
    size_t t=maxVertex+2;

    verticesInGroup=std::vector<size_t>();
    verticesInGroup.push_back(s);
    vToGroup.push_back(0);
    groups[0]=verticesInGroup;

    verticesInGroup=std::vector<size_t>();
    verticesInGroup.push_back(t);
    vToGroup.push_back(frameCounter);
    groups[frameCounter]=verticesInGroup;


    //TODO first and last frame for s and t!

}


template<class T>
template<class STR>
inline void VertexGroups<T>::initFromStream(STR& timeData,size_t maxTimeToRead){
	size_t lineCounter=0;
	std::vector<size_t> currentGroup;
	std::vector<std::string> strings;
	std::string line;
	char delim=',';

	unsigned int previousTime=0;
	unsigned int time=0;

	while (std::getline(timeData, line) && !line.empty()) {

		//std::cout<<"time file line "<<lineCounter<<" delim "<<delim<<std::endl;
		lineCounter++;

		strings = split(line,delim);
		//std::cout<<"split "<<std::endl;
		if (strings.size() < 2) {
			throw std::runtime_error(
					std::string("Vertex and time frame expected"));
		}

		unsigned int v = std::stoul(strings[0]);
		time = std::stoul(strings[1]);
		if(time>maxTimeToRead){
			//groups[previousTime]=currentGroup;
			break;
		}
		if(vToGroup.size()!=v){
			throw std::runtime_error(
					std::string("Wrong vertex numbering in time file"));
		}
		else{

			vToGroup.push_back(time);
			//if(time==groups.size()+1){  //Change groups to map?
			if(time==previousTime||previousTime==0){
				currentGroup.push_back(v);
			}
			//else if(time==groups.size()+2){
			else{
				//groups.push_back(currentGroup);
				groups[previousTime]=currentGroup;
				currentGroup=std::vector<size_t>();
				currentGroup.push_back(v);
			}
			previousTime=time;


		}
	}
	groups[previousTime]=currentGroup;


	maxTime=*(vToGroup.rbegin());

	//time frame of s is zero
	vToGroup.push_back(0);
	currentGroup=std::vector<size_t>();
	currentGroup.push_back(vToGroup.size()-1);
	groups[0]=currentGroup;

	//time frame of t is maxTime
	vToGroup.push_back(maxTime+1);
	currentGroup=std::vector<size_t>();
	currentGroup.push_back(vToGroup.size()-1);
	groups[maxTime+1]=currentGroup;

	maxVertex=vToGroup.size()-3;
}

//
//template<class T>
//std::vector<std::vector<size_t>> extractInnerPaths(VertexGroups<size_t>& vg,std::vector<std::vector<size_t>>& paths,T minT,T maxT){
//	//maxT inclusive
//	std::vector<std::vector<size_t>> outputPaths;
//	for (int i = 0; i < paths.size(); ++i) {
//		std::vector<size_t>& path=paths[i];
//		size_t counter=0;
//		bool startFound=false;
//		while(!startFound&&counter<path.size()){
//			size_t time=vg.getGroupIndex(path[counter]);
//			if(time>=minT){
//				startFound=true;
//			}
//			else{
//				counter++;
//			}
//		}
//		if(startFound){
//			std::vector<size_t> newPath;
//			bool endFound=false;
//			while(!endFound&&counter<path.size()){
//
//				size_t time=vg.getGroupIndex(path[counter]);
//				if(time<=maxT){
//					newPath.push_back(path[counter]);
//					counter++;
//				}
//				else{
//					endFound=true;
//				}
//			}
//			outputPaths.push_back(newPath);
//		}
//
//	}
//	return outputPaths;
//}




template<class T>
std::vector<std::vector<size_t>> extractInnerPaths(VertexGroups<size_t>& vg,std::vector<std::vector<size_t>>& paths,T minT,T maxT){
	//maxT inclusive
	std::vector<std::vector<size_t>> outputPaths;
	for (int i = 0; i < paths.size(); ++i) {
		std::vector<size_t>& path=paths[i];
		std::vector<size_t> outputPath;
		for(size_t vertex: path){
			size_t time=vg.getGroupIndex(vertex);
			if(time<=maxT){
				if(time>=minT){
					outputPath.push_back(vertex);
				}
			}
			else break;
		}
		if(outputPath.size()>0){
			outputPaths.push_back(outputPath);
		}
	}
	return outputPaths;
}




template<class T>
	std::vector<std::unordered_set<size_t>> initReachableSet(T & graph,DisjointParams<>& parameters,VertexGroups<size_t>* vg=0){

	levinkov::Timer tfw;
			tfw.start();
		std::cout<<"Run Floyd Warshall"<<std::endl;
		parameters.infoFile()<<"Run Floyd Warshall"<<std::endl;
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

		std::cout << "fw finished in time "<<tfw.get_elapsed_seconds() << std::endl;
		parameters.infoFile()<<"fw finished in time "<<tfw.get_elapsed_seconds() << std::endl;
		return desc;

}



template<class T>
	std::vector<std::vector<bool>> initReachable(T & graph,DisjointParams<>& parameters,VertexGroups<size_t>* vg=0){

	levinkov::Timer tfw;
			tfw.start();
		std::cout<<"Run Floyd Warshall"<<std::endl;
		parameters.infoFile()<<"Run Floyd Warshall"<<std::endl;
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

		std::cout << "fw finished in time "<<tfw.get_elapsed_seconds() << std::endl;
		parameters.infoFile()<<"fw finished in time "<<tfw.get_elapsed_seconds() << std::endl;
		return desc;

}





template<class T = size_t>
struct CompleteStructure {
public:
	CompleteStructure(DisjointParams<T> & configParameters,char delim=','):
		params(configParameters)
{
		vg=VertexGroups<>(configParameters,delim);
		maxTime=vg.getMaxTime();
		std::cout<<"max time "<<maxTime<<std::endl;
		configParameters.infoFile()<<"max time "<<maxTime<<std::endl;
		completeGraph=andres::graph::Digraph<>(vg.getMaxVertex()+1);
		std::cout<<"cg vertices "<<completeGraph.numberOfVertices()<<std::endl;
		configParameters.infoFile()<<"cg vertices "<<completeGraph.numberOfVertices()<<std::endl;

}

	void addEdgesFromFile();
    void addEdges(size_t time1,size_t time2,const py::array_t<double> inputMatrix);
    //void addEdges(size_t time1,size_t time2,std::stringstream& data);

	andres::graph::Digraph<> completeGraph;
	std::vector<double> completeScore;
	VertexGroups<> vg;
	size_t maxTime;
	DisjointParams<T>& params;

};

template<class T>
inline void CompleteStructure<T>::addEdges(size_t time1,size_t time2,const py::array_t<double> inputMatrix){

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
            size_t vertex2=j+transformIndex2;
            completeGraph.insertEdge(vertex1,vertex2);
            completeScore.push_back(score);
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
inline void CompleteStructure<T>::addEdgesFromFile(){
	char delim=',';

	std::string line;
	std::ifstream data;
	try{
		data.open(params.getGraphFileName());
		if(!data){
			throw std::system_error(errno, std::system_category(), "failed to open "+params.getGraphFileName());
		}

		std::getline(data, line);
		double objValue=0;

		std::cout << "Read big graph" << std::endl;
		params.infoFile()<<"Read big graph" << std::endl;
		std::vector<std::string> strings;


		std::cout<<"Reading vertices from file. "<<std::endl;
		params.infoFile()<<"Reading vertices from file. "<<std::endl;
		params.infoFile().flush();
		//Vertices that are not found have score=0. Appearance and disappearance cost are read here.
		while (std::getline(data, line) && !line.empty()) {

		}

		std::cout<<"Reading base edges from file. "<<std::endl;
		params.infoFile()<<"Reading base edges from file. "<<std::endl;
		params.infoFile().flush();
		size_t maxLabel=0;
		while (std::getline(data, line) && !line.empty()) {

			strings = split(line, delim);

			unsigned int v = std::stoul(strings[0]);
			unsigned int w = std::stoul(strings[1]);
			size_t l0=vg.getGroupIndex(v);
			size_t l1=vg.getGroupIndex(w);
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
}


#endif /* INCLUDE_DISJOINT_PATHS_DISJOINTPATHSMETHODS_HXX_ */
