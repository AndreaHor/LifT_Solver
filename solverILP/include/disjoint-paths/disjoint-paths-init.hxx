/*
 * disjoint-paths-init.hxx
 *
 *  Created on: Sep 10, 2018
 *      Author: fuksova
 */

#ifndef INCLUDE_DISJOINT_PATHS_DISJOINT_PATHS_INIT_HXX_
#define INCLUDE_DISJOINT_PATHS_DISJOINT_PATHS_INIT_HXX_

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
#include "disjoint-paths/inputConfig.hxx"
#include <utility>

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

	VertexGroups(disjointPaths::ConfigDisjoint<>& parameters,char delim){
		//std::vector<std::vector<size_t>> groups;

		std::string line;
		std::vector<std::string> strings;
		std::ifstream timeData;
		try{
			timeData.open(parameters.getTimeFileName());
			if(!timeData){
				throw std::system_error(errno, std::system_category(), "failed to open "+parameters.getTimeFileName());
			}

			size_t lineCounter=0;
			std::vector<size_t> currentGroup;

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
				if(time>parameters.getMaxTimeFrame()){
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
					//				else{
					//					throw std::runtime_error("Wrong time frame numbering.");
					//				}

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

			//		for (int i = 0; i < maxTime+2; ++i) {
			//			std::cout<<"group "<<i<<": ";
			//			for (int j = 0; j < groups[i].size(); ++j) {
			//				std::cout<<groups[i][j]<<",";
			//			}
			//			std::cout<<std::endl;
			//
			//		}
		}
		catch (std::system_error& er) {
			std::clog << er.what() << " (" << er.code() << ")" << std::endl;

		}

	}



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
	std::vector<std::unordered_set<size_t>> initReachableSet(T & graph,ConfigDisjoint<>& parameters,VertexGroups<size_t>* vg=0){

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
	std::vector<std::vector<bool>> initReachable(T & graph,ConfigDisjoint<>& parameters,VertexGroups<size_t>* vg=0){

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

	CompleteStructure(ConfigDisjoint<T> & configParameters,char delim=','){

		vg=VertexGroups<>(configParameters,delim);
		maxTime=vg.getMaxTime();
		std::cout<<"max time "<<maxTime<<std::endl;
		configParameters.infoFile()<<"max time "<<maxTime<<std::endl;
		completeGraph=andres::graph::Digraph<>(vg.getMaxVertex()+1);
		std::cout<<"cg vertices "<<completeGraph.numberOfVertices()<<std::endl;
		configParameters.infoFile()<<"cg vertices "<<completeGraph.numberOfVertices()<<std::endl;
		std::string line;
		std::ifstream data;
		try{
			data.open(configParameters.getGraphFileName());
			if(!data){
				throw std::system_error(errno, std::system_category(), "failed to open "+configParameters.getGraphFileName());
			}

			std::getline(data, line);
			double objValue=0;

			std::cout << "Read big graph" << std::endl;
			configParameters.infoFile()<<"Read big graph" << std::endl;
			std::vector<std::string> strings;


			std::cout<<"Reading vertices from file. "<<std::endl;
			configParameters.infoFile()<<"Reading vertices from file. "<<std::endl;
			configParameters.infoFile().flush();
			//Vertices that are not found have score=0. Appearance and disappearance cost are read here.
			while (std::getline(data, line) && !line.empty()) {

			}

			std::cout<<"Reading base edges from file. "<<std::endl;
			configParameters.infoFile()<<"Reading base edges from file. "<<std::endl;
			configParameters.infoFile().flush();
			size_t maxLabel=0;
			while (std::getline(data, line) && !line.empty()) {

				strings = split(line, delim);

				unsigned int v = std::stoul(strings[0]);
				unsigned int w = std::stoul(strings[1]);
				size_t l0=vg.getGroupIndex(v);
				size_t l1=vg.getGroupIndex(w);
				if(v>vg.getMaxVertex()||w>vg.getMaxVertex()) continue;
				//std::cout<<"edge "<<v<<" "<<w<<std::endl;

				if(l1-l0<=configParameters.getMaxTimeGapComplete()){
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

	andres::graph::Digraph<> completeGraph;
	std::vector<double> completeScore;
	VertexGroups<> vg;
	size_t maxTime;

};



template<class T = size_t>
struct DisjointStructure {
public:
	typedef T size_type;


	DisjointStructure(disjointPaths::ConfigDisjoint<>& configParameters,char delim=',',CompleteStructure<>* cs=0,size_t minTime=0,size_t maxTime=0);
	bool isReachable(size_t i,size_t j) const{
		if(i==t_) return false;  //Assume no path from the terminal node
		if(desc.size()==0) return true;
		return desc[i][j];


	}

	bool isReachableNew(size_t i,size_t j) const{
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


	std::vector<std::vector<bool>>* getPReachable(){
		return &desc;
	}

	std::vector<std::unordered_set<size_t>>* getPReachableNew(){
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

	const std::vector<double>& getLiftedEdgesScore() {
			return liftedEdgeScore;
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


	disjointPaths::ConfigDisjoint<>& parameters;

	bool isTimeFramesEnabled(){
		return useTimeFrames;
	}

	VertexGroups<size_t> vertexGroups;
	size_t minV=0;
	size_t maxV=0;


private:

	size_t s_;
	size_t t_;

	std::vector<double> vertexScore;
	std::vector<double> edgeScore;
	std::vector<double> liftedEdgeScore;
	std::vector<std::vector<bool>> desc;
	std::vector<std::unordered_set<size_t>> reachable;



	andres::graph::Digraph<> graph_;
	andres::graph::Digraph<> graphLifted_;


	void readGraph(std::ifstream& data,size_t maxVertex,char delim);
	void readGraphWithTime(size_t minTime,size_t maxTime,CompleteStructure<>* cs);
	//void persistentEdges();
	void createKnnBaseGraph();
	void keepFractionOfLiftedNew();
	void keepFractionOfLifted();
	void keepFractionOfLifted3();
	bool useTimeFrames;


};


template<class T>
inline DisjointStructure<T>::DisjointStructure(disjointPaths::ConfigDisjoint<>& configParameters,char delim,CompleteStructure<>* cs,size_t minTime,size_t maxTime):
parameters(configParameters)
{

	std::cout<<"interval "<<minTime<<","<<maxTime<<std::endl;
	parameters.infoFile()<<"interval "<<minTime<<","<<maxTime<<std::endl;
	parameters.infoFile().flush();

	useTimeFrames=parameters.isRestrictFrames()||parameters.isSparsify();
	size_t maxVertex;
	if(cs==0){
	//	if(useTimeFrames){
			vertexGroups=VertexGroups<size_t>(configParameters,delim);
			maxVertex=vertexGroups.getMaxVertex();
//		}
//		else{
//			maxVertex=std::numeric_limits<size_t>::max();
//		}
		std::ifstream graphFile;
		try{
			graphFile.open(parameters.getGraphFileName());
			if(!graphFile){
				throw std::system_error(errno, std::system_category(), "failed to open "+parameters.getGraphFileName());
			}
			readGraph(graphFile,maxVertex,delim);
			if(!parameters.isAutomaticLifted()){
				//std::vector<std::vector<bool>> secOrderDesc=automaticLifted(graph_);
				std::cout<<"Reading lifted edges from file."<<std::endl;
				std::string line;
				std::vector<std::string> strings;
				while (std::getline(graphFile, line) && !line.empty()) {
					strings = split(line, delim);
					if (strings.size() < 3) {
						throw std::runtime_error(
								std::string("Edge vertices and score expected on every edge line "));
					}

					unsigned int v = std::stoul(strings[0]);
					unsigned int w = std::stoul(strings[1]);
					if(v>=graph_.numberOfVertices()-2||w>=graph_.numberOfVertices()-2) continue;
					//if(isReachable(v,w)&&v!=s_&&w!=t_&&v!=t_&&w!=s_){
					if(isReachableNew(v,w)&&v!=s_&&w!=t_&&v!=t_&&w!=s_){
						double score = std::stod(strings[2]);
						//if(secOrderDesc[v][w]){
						auto edgeTest=graphLifted_.findEdge(v,w);
						if(!edgeTest.first){
							graphLifted_.insertEdge(v, w);
							liftedEdgeScore.push_back(score);
						}
						else{
							liftedEdgeScore[edgeTest.second]=score;
						}
					}
				}

			}
			graphFile.close();
		}
		catch (std::system_error& er) {
			std::clog << er.what() << " (" << er.code() << ")" << std::endl;

		}
	}
	else{
		readGraphWithTime(minTime,maxTime,cs);

	}


	//desc=disjointPaths::initReachable(graph_);

	//std::vector<std::vector<bool>> secOrderDesc=automaticLifted(graph_);

	//TODO ensure that the lifted edges are added only if the reachability and non uniqueness is satisfied!

	if(graph_.numberOfVertices()>2){

		if(parameters.getMaxTimeLifted()>0&&parameters.isAutomaticLifted()){
			std::cout<<"Adding automatic lifted edges"<<std::endl;
			parameters.infoFile()<<"Adding automatic lifted edges"<<std::endl;
			parameters.infoFile().flush();
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
			std::cout<<"done"<<std::endl;
			parameters.infoFile()<<"done"<<std::endl;
			parameters.infoFile().flush();

		}


		std::cout<<"number of vertices "<<graph_.numberOfVertices()<<std::endl;
		parameters.infoFile()<<"number of vertices "<<graph_.numberOfVertices()<<std::endl;
		parameters.infoFile().flush();
		if(parameters.isSparsify()){
			createKnnBaseGraph();

			if(parameters.getMaxTimeLifted()>0){
				if(parameters.isNewLifted()){
					keepFractionOfLifted3();
				}
				else{
					keepFractionOfLifted();
				}
			}
		}
		else{
			//desc=initReachable(graph_,parameters);
			reachable=initReachableSet(graph_,parameters);
		}
	}
}


template<class T>
inline void DisjointStructure<T>::readGraphWithTime(size_t minTime,size_t maxTime,CompleteStructure<>* cs){

	andres::graph::Digraph<>& completeGraph=cs->completeGraph;
	std::vector<double>& completeScore=cs->completeScore;
	VertexGroups<> vg=cs->vg;

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
		std::cout<<"empty interval in ds constructor"<<std::endl;

	}
	else{
		std::cout<<"valid interval in ds constructor"<<std::endl;
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
		std::cout<<"DS constructor max Vertex "<<maxV<<std::endl;
		//}
	}

}







template<class T>
inline void DisjointStructure<T>::readGraph(std::ifstream& data,size_t maxVertex,char delim){
	std::string line;
	//	char delim = ' ';
	size_t lineCounter=0;
	std::getline(data, line);
	lineCounter++;
	std::cout << "called read graph" << std::endl;
	parameters.infoFile()<<"called read graph" << std::endl;
	std::vector<std::string> strings = split(line, delim);
	size_t numberOfVertices;

	if (strings.size() == 1) {
		if(parameters.isRestrictFrames()){
			numberOfVertices=maxVertex+3;
		}
		else{
			numberOfVertices = stoul(strings[0]);
			numberOfVertices += 2;  //to include s and t in the end of the list
		}
		s_ = numberOfVertices - 2;
		t_ = numberOfVertices - 1;

	} else {
		std::string str="first row must contain 1 number, detected ";
		str+=std::to_string(strings.size());
		str+="numbers";
		throw std::runtime_error(str);
	}



	graphLifted_ = andres::graph::Digraph<>(numberOfVertices);
	graph_ = andres::graph::Digraph<>(numberOfVertices);
	std::vector<double> inputCosts(numberOfVertices-2,parameters.getInputCost());
	std::vector<double> outputCosts(numberOfVertices-2,parameters.getOutputCost());


	// std::vector<std::pair<size_t,siz  Data<>e_t> > liftedEdges;
	vertexScore = std::vector<double>(numberOfVertices, 0);

	std::cout<<"Reading vertices from file. "<<std::endl;
	parameters.infoFile()<<"Reading vertices from file. "<<std::endl;
	parameters.infoFile().flush();
	//Vertices that are not found have score=0. Appearance and disappearance cost are read here.
	bool vertexAndInOutScore=false;
	while (std::getline(data, line) && !line.empty()) {
		if(vertexAndInOutScore){ //By default disabled in the final version
			lineCounter++;
			strings = split(line, delim);
			if (strings.size() < 2) {
				throw std::runtime_error(
						std::string("Vertex and its score expected"));
			}


			unsigned int v = std::stoul(strings[0]);
			if(v>graph_.numberOfVertices()-3) continue;
			double score = std::stod(strings[1]);
			vertexScore[v] = score;

			if(strings.size()==4){
				inputCosts[v]=std::stod(strings[2]);
				outputCosts[v]=std::stod(strings[3]);
			}
		}

	}

	for (int v = 0; v < numberOfVertices-2; ++v) {
		graph_.insertEdge(s_,v);
		edgeScore.push_back(inputCosts[v]);
		graph_.insertEdge(v,t_);
		edgeScore.push_back(outputCosts[v]);
	}

	size_t maxGap=parameters.getMaxTimeGapComplete();

	std::cout<<"Reading base edges from file. "<<std::endl;
	parameters.infoFile()<<"Reading base edges from file. "<<std::endl;
	parameters.infoFile().flush();
	while (std::getline(data, line) && !line.empty()) {
		lineCounter++;
		strings = split(line, delim);
		if (strings.size() < 3) {
			throw std::runtime_error(
					std::string("Edge vertices and score expected, line "+std::to_string(lineCounter)));
		}

		unsigned int v = std::stoul(strings[0]);
		unsigned int w = std::stoul(strings[1]);

		if(v>numberOfVertices-3||w>numberOfVertices-3) continue;

		size_t gv=vertexGroups.getGroupIndex(v);
		size_t gw=vertexGroups.getGroupIndex(w);
		if(gw-gv>maxGap) continue;

		//if(v>=graph_.numberOfVertices()-2||w>=graph_.numberOfVertices()-2) continue;
		double score = std::stod(strings[2]);
		auto edgeTest=graph_.findEdge(v,w);

		if(!edgeTest.first){  //if the edge does not exist
			graph_.insertEdge(v, w);
			edgeScore.push_back(score);

		}
		else{  //if the edge already exists, only update the score
			edgeScore[edgeTest.second]=score;

		}


	}

}







template<class T>
inline void DisjointStructure<T>::createKnnBaseGraph(){
	std::cout<<"Sparsify base graph"<<std::endl;
	parameters.infoFile()<<"Sparsify base graph"<<std::endl;
	andres::graph::Digraph<> tempGraph(graph_.numberOfVertices());
	std::vector<double> newBaseCosts;
	//std::vector<size_t> inOutEdges;
	size_t k=parameters.getKnnK();
	//std::vector<size_t> goodLongEdges;


	std::vector<bool> finalEdges(graph_.numberOfEdges(),false);
	for (int v0 = 0; v0 < graph_.numberOfVertices(); ++v0) {
		std::unordered_map<int,std::list<size_t>> edgesToKeep;
		size_t l0=vertexGroups.getGroupIndex(v0);
		for (size_t ne = 0; ne < graph_.numberOfEdgesFromVertex(v0); ++ne) {
			size_t e=graph_.edgeFromVertex(v0,ne);
			size_t v1=graph_.vertexFromVertex(v0,ne);
//			std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<problemGraph.getEdgeScore(e)<<std::endl;
			if(v0==s_||v1==t_){
				//tempGraph.insertEdge(v0,v1);
				//newBaseCosts.push_back(edgeScore[e]);
				finalEdges[e]=true;
			}
			else{
				size_t l1=vertexGroups.getGroupIndex(v1);
				size_t gap=l1-l0;
				if(gap<=parameters.getMaxTimeBase()){
				//if(gap<=parameters.getKnnTimeGap()){
					//gap=std::min(parameters.getKnnTimeGap()+1,gap);
					double cost=edgeScore[e];
					if(edgesToKeep.count(gap)>0){
						std::list<size_t>& smallList=edgesToKeep[gap];
						auto it=smallList.begin();
						double bsf=edgeScore[*it];
						//std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<bsf<<std::endl;
						while(bsf>cost&&it!=smallList.end()){
							it++;
							size_t index=*it;
							if(it!=smallList.end()){
								bsf=edgeScore[index];
								//	std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<bsf<<std::endl;
							}
						}
						if(it!=smallList.begin()){
							smallList.insert(it,e);
							if(smallList.size()>k) smallList.pop_front();
						}
						else if(smallList.size()<k){
							smallList.push_front(e);
						}
					}
					else{
						edgesToKeep[gap].push_front(e);
					}
				}
//				else if(gap<=parameters.getMaxTimeBase()){
//					if(getEdgeScore(e)<=parameters.getBaseUpperThreshold()){
//						//tempGraph.insertEdge(v0,v1);
//						//newBaseCosts.push_back(getEdgeScore(e));
//						finalEdges[e]=true;
//					}
//
//				}
			}
		}
		//std::cout.precision(4);
		double bsf=0;
		for (int gap = 0; gap <= parameters.getKnnTimeGap(); ++gap) {
			if(edgesToKeep.count(gap)>0){
				auto& smallList=edgesToKeep[gap];
				for(size_t e:smallList){
					finalEdges[e]=true;
					if(edgeScore[e]<bsf){
						bsf=edgeScore[e];
					}
				}
			}
		}

		bool onlyImproving=parameters.isRequireImproving();
		for (int gap =  parameters.getKnnTimeGap()+1;gap<=parameters.getMaxTimeBase(); ++gap) {
			if(edgesToKeep.count(gap)>0){
				if(onlyImproving){
					double currentBsf=0;
					auto& smallList=edgesToKeep[gap];
					for(size_t e:smallList){
						double score=edgeScore[e];
						if(score<bsf&&score<=parameters.getBaseUpperThreshold()){
							finalEdges[e]=true;
							if(score<currentBsf){
								currentBsf=score;
							}
						}
					}
					bsf=std::min(bsf,currentBsf);
				}
				else{
					auto& smallList=edgesToKeep[gap];
					for(size_t e:smallList){
						double score=edgeScore[e];
						if(score<=parameters.getBaseUpperThreshold()){
							finalEdges[e]=true;
						}
					}
				}
			}
		}

//		for (auto itMap=edgesToKeep.begin();itMap!=edgesToKeep.end();itMap++) {
//			std::list<size_t>& smallList=itMap->second;
//			//for(std::list<size_t>::reverse_iterator itList=smallList.rbegin();itList!=smallList.rend();++itList){
//			for(auto itList=smallList.begin();itList!=smallList.end();++itList){
//				size_t edgeIndex=*itList;
//				finalEdges[edgeIndex]=true;
//
//
//			}
//		}
	}

	for (int e = 0; e < graph_.numberOfEdges(); ++e) {
		if(finalEdges[e]){
			size_t v0=graph_.vertexOfEdge(e,0);
			size_t v1=graph_.vertexOfEdge(e,1);
			tempGraph.insertEdge(v0,v1);
			newBaseCosts.push_back(edgeScore[e]);
		}
	}

	if(newBaseCosts.size()!=tempGraph.numberOfEdges()){
		throw std::runtime_error("Error in base graph sparsification.");
	}



	graph_=tempGraph;
	edgeScore=newBaseCosts;

	if(graph_.numberOfEdges()!=newBaseCosts.size()){
		parameters.infoFile()<<"edge number mismatch, graph: "<<graph_.numberOfEdges()<<", cost vector "<<newBaseCosts.size()<<std::endl;
		parameters.infoFile().flush();
		std::cout<<"edge number mismatch, graph: "<<graph_.numberOfEdges()<<", cost vector "<<newBaseCosts.size()<<std::endl;
	}
	else{
		std::cout<<"edge number and graph size match "<<std::endl;
		parameters.infoFile()<<"edge number and graph size match "<<std::endl;
		parameters.infoFile().flush();
	}


	if(parameters.getMaxTimeLifted()>0){
		if(parameters.getSmallIntervals()==0){
			//desc=initReachable(graph_,parameters,&vertexGroups);
			reachable=initReachableSet(graph_,parameters,&vertexGroups);
		}
		//desc=initReachable(graph_,parameters);
		reachable=initReachableSet(graph_,parameters);


		std::cout<<"Left "<<newBaseCosts.size()<<" base edges"<<std::endl;
		parameters.infoFile()<<"Left "<<newBaseCosts.size()<<" base edges"<<std::endl;
		parameters.infoFile().flush();
	}
	else{
		std::cout<<"FW skipped"<<std::endl;
	}

}


template<class T>
inline void DisjointStructure<T>::keepFractionOfLifted(){

		std::cout<<"Sparsify lifted graph"<<std::endl;
		parameters.infoFile()<<"Sparsify lifted graph"<<std::endl;
		parameters.infoFile().flush();
		//TODO run automaticLifted to find candidates first

		double negMaxValue=0;
		double posMinValue=0;
		bool useAdaptive=false;
		//TODO adaptive lifted threshold in config file
		std::vector<double> newLiftedCosts;

		andres::graph::Digraph<> tempGraphLifted=(graph_.numberOfVertices());


		if(useAdaptive){

			double fPositive=parameters.getPositiveThresholdLifted();
			double fNegative=parameters.getNegativeThresholdLifted();

			assert(fNegative>=0&&fNegative<=1);
			assert(fPositive>=0&&fPositive<=1);

			std::priority_queue<double> allPositive;
			std::priority_queue<double> allNegative;
			for (int i = 0; i < graphLifted_.numberOfEdges(); ++i) {
				double cost=getLiftedEdgeScore(i);
				if(cost<0){
					allNegative.push(-cost);
				}
				else{
					allPositive.push(cost);
				}
			}

			double negativeToKeep=fNegative*allNegative.size();
			double positiveToKeep=fPositive*allPositive.size();
			int negNumber=int(round(negativeToKeep));
			int posNumber=int(round(positiveToKeep));

			for (int i = 0; i < negNumber; ++i) {
				allNegative.pop();
			}
			for (int i = 0; i < posNumber; ++i) {
				allPositive.pop();
			}



			if(!allNegative.empty()){
				negMaxValue=-allNegative.top();
			}
			std::cout<<"max negative lifted "<<negMaxValue<<std::endl;
			parameters.infoFile()<<"max negative lifted "<<negMaxValue<<std::endl;

			if(!allPositive.empty()){
				posMinValue=allPositive.top();
			}
			std::cout<<"min positive lifted "<<posMinValue<<std::endl;
			parameters.infoFile()<<"min positive lifted "<<posMinValue<<std::endl;

		}
		else{
			negMaxValue=parameters.getNegativeThresholdLifted();
			posMinValue=parameters.getPositiveThresholdLifted();
		}

//		std::unordered_map<size_t,std::set<size_t>> liftedEdges;
//			for (int v = 0; v < graphLifted_.numberOfVertices()-2; ++v) {
//				std::vector<bool> alternativePath(graphLifted_.numberOfVertices(),false);
//				for (int i = 0; i < graph_.numberOfEdgesFromVertex(v); ++i) {
//					size_t w=graph_.vertexFromVertex(v,i);
//					for (int u = 0; u < graph_.numberOfVertices(); ++u) {
//						//if(isReachable(w,u)&&w!=u) alternativePath[u]=true;
//						if(isReachableNew(w,u)&&w!=u) alternativePath[u]=true;
//					}
//				}
//				for (int i = 0; i < graphLifted_.numberOfEdgesFromVertex(v); ++i) {
//					size_t w=graphLifted_.vertexFromVertex(v,i);
//					if(w!=t_){
//						if(alternativePath[w]) liftedEdges[v].insert(w);
//
//					}
//
//				}
//			}


		std::unordered_map<size_t,std::set<size_t>> liftedEdges;
		for (int v = 0; v < graphLifted_.numberOfVertices()-2; ++v) {
			std::unordered_set<size_t> alternativePath;
			for (int i = 0; i < graph_.numberOfEdgesFromVertex(v); ++i) {
				size_t w=graph_.vertexFromVertex(v,i);
				for(size_t u:reachable[w]){
					if(u!=w) alternativePath.insert(u);
				}
			}
			for (int i = 0; i < graphLifted_.numberOfEdgesFromVertex(v); ++i) {
				size_t w=graphLifted_.vertexFromVertex(v,i);
				if(w!=t_){
					if(alternativePath.count(w)>0) liftedEdges[v].insert(w);

				}

			}
		}


			std::cout<<"done"<<std::endl;
			parameters.infoFile()<<"done"<<std::endl;
			parameters.infoFile().flush();


		for (int i = 0; i < graphLifted_.numberOfEdges(); ++i) {
			size_t v0=graphLifted_.vertexOfEdge(i,0);
			size_t v1=graphLifted_.vertexOfEdge(i,1);
			int l0=vertexGroups.getGroupIndex(v0);
			int l1=vertexGroups.getGroupIndex(v1);
			double cost=getLiftedEdgeScore(i);
			bool goodCost=(cost<negMaxValue)||(cost>posMinValue);
			//if(isReachable(v0,v1)){
			if(isReachableNew(v0,v1)){

				int timeGapDiff=l1-l0-parameters.getDenseTimeLifted();
				bool timeConstraint=l1-l0<=parameters.getDenseTimeLifted()||((l1-l0)<=parameters.getMaxTimeLifted()&&(timeGapDiff%parameters.getLongerIntervalLifted())==0);
				//int modulo=timeGapDiff%parameters.getLongerIntervalLifted();
				//if(l1-l0<=parameters.getDenseTimeLifted()||((l1-l0)<=parameters.getMaxTimeLifted()&&(timeGapDiff%parameters.getLongerIntervalLifted())==0)){
				if(timeConstraint&&goodCost){
					if(liftedEdges[v0].count(v1)>0){
						tempGraphLifted.insertEdge(v0,v1);
						newLiftedCosts.push_back(cost);
					}
					else{
						auto edgeTest=graph_.findEdge(v0,v1);
						if(edgeTest.first){
							edgeScore[edgeTest.second]+=cost;  //Compensate that the lifted edge has been removed
						}
//						else{
//							throw std::runtime_error("Error in lifted graph sparsification.");
//						}
					}

				}
			}

		}






		liftedEdgeScore=newLiftedCosts;

		graphLifted_=tempGraphLifted;
		std::cout<<"Left "<<newLiftedCosts.size()<<" lifted edges."<<std::endl;
		parameters.infoFile()<<"Left "<<newLiftedCosts.size()<<" lifted edges."<<std::endl;
		parameters.infoFile().flush();

		if(graphLifted_.numberOfEdges()!=newLiftedCosts.size()){
			std::cout<<"lifted edge number mismatch, lifted graph: "<<graphLifted_.numberOfEdges()<<", cost vector "<<newLiftedCosts.size()<<std::endl;
			parameters.infoFile()<<"lifted edge number mismatch, lifted graph: "<<graphLifted_.numberOfEdges()<<", cost vector "<<newLiftedCosts.size()<<std::endl;
		}
		else{
			std::cout<<"lifted edge number and lifted graph size match "<<std::endl;
			parameters.infoFile()<<"lifted edge number and lifted graph size match "<<std::endl;

		}
		parameters.infoFile().flush();

	}





template<class T>
inline void DisjointStructure<T>::keepFractionOfLifted3(){
	std::cout<<"Sparsify lifted graph"<<std::endl;
	parameters.infoFile()<<"Sparsify lifted graph"<<std::endl;
	parameters.infoFile().flush();

	double negMaxValue=parameters.getNegativeThresholdLifted();
	double posMinValue=parameters.getPositiveThresholdLifted();

	size_t interval=parameters.getLongerIntervalLifted();

	size_t numberOfSparseGaps=parameters.getMaxTimeLifted()-parameters.getDenseTimeLifted();
	size_t numberOfIntervals=(numberOfSparseGaps)/interval;
	if((numberOfSparseGaps)%interval!=0){
		numberOfIntervals++;
	}

	//TODO make this more effective
	std::unordered_map<size_t,std::set<size_t>> liftedEdges;
	for (int v = 0; v < graphLifted_.numberOfVertices()-2; ++v) {
		std::vector<bool> alternativePath(graphLifted_.numberOfVertices(),false);
		for (int i = 0; i < graph_.numberOfEdgesFromVertex(v); ++i) {
			size_t w=graph_.vertexFromVertex(v,i);
			for (int u = 0; u < graph_.numberOfVertices(); ++u) {
				//if(isReachable(w,u)&&w!=u) alternativePath[u]=true;
				if(isReachableNew(w,u)&&w!=u) alternativePath[u]=true;
			}
		}
		for (int i = 0; i < graphLifted_.numberOfEdgesFromVertex(v); ++i) {
			size_t w=graphLifted_.vertexFromVertex(v,i);
			if(w!=t_&&alternativePath[w]) liftedEdges[v].insert(w);
		}
	}
	std::cout<<"second order reachability done"<<std::endl;
	parameters.infoFile()<<"second order reachability done"<<std::endl;
	parameters.infoFile().flush();



	std::vector<bool> finalEdges(graphLifted_.numberOfEdges(),false);

	for (int v = 0; v < graphLifted_.numberOfVertices()-2; ++v) {
		std::map<size_t,std::vector<size_t>> gapToPositiveEdges;
		std::map<size_t,std::vector<size_t>> gapToNegativeEdges;
		std::map<size_t,double> gapToPosCost;
		std::map<size_t,double> gapToNegCost;
        int l0=vertexGroups.getGroupIndex(v);

        for (int i = 0; i < graphLifted_.numberOfEdgesFromVertex(v); ++i) {
        	size_t edge=graphLifted_.edgeFromVertex(v,i);
        	size_t w=graphLifted_.vertexFromVertex(v,i);

        	int l1=vertexGroups.getGroupIndex(w);
        	int gap=l1-l0;
        	//if(isReachable(v,w)&&gap<=parameters.getMaxTimeLifted()){
        	if(isReachableNew(v,w)&&gap<=parameters.getMaxTimeLifted()){
        		double cost=getLiftedEdgeScore(edge);
        		if(gap<=parameters.getDenseTimeLifted()){
        			if((cost<=negMaxValue)||(cost>=posMinValue)){
        				finalEdges[edge]=true;
        			}
        		}
        		else{
        			if(cost<=negMaxValue){
        				gapToNegativeEdges[gap].push_back(edge);
        				if(gapToNegCost.count(gap)>1){
        					gapToNegCost[gap]+=cost;
        				}
        				else{
        					gapToNegCost[gap]=cost;
        				}

        			}
        			else if(cost >=posMinValue){
        				gapToPositiveEdges[gap].push_back(edge);
        				if(gapToPosCost.count(gap)>0){
        					gapToPosCost[gap]+=cost;
        				}
        				else{
        					gapToPosCost[gap]=cost;
        				}

        			}
        		}
        	}


        }


        for (int i = 0; i < numberOfIntervals; ++i) {
        	double positiveBSF=0;
        	double negativeBSF=0;
        	size_t indexPositiveBSF=i*interval+1;
        	size_t indexNegativeBSF=i*interval+1;

        	for (int j = 1; j <= parameters.getLongerIntervalLifted(); ++j) {
        		size_t gap=i*interval+j;
        		if(gap>parameters.getMaxTimeLifted()) break;
        		if(positiveBSF<gapToPosCost[gap]){
        			positiveBSF=gapToPosCost[gap];
        			indexPositiveBSF=gap;
        		}
        		if(negativeBSF>gapToNegCost[gap]){
        			negativeBSF=gapToNegCost[gap];
        			indexNegativeBSF=gap;
        		}
        	}
        	for(size_t e:gapToPositiveEdges[indexPositiveBSF]){
        		finalEdges[e]=true;
        	}
        	for(size_t e:gapToNegativeEdges[indexNegativeBSF]){
        		finalEdges[e]=true;
        	}
        }

	}

	std::vector<double> newLiftedCosts;
	andres::graph::Digraph<> tempGraphLifted=(graph_.numberOfVertices());

	for (int e = 0; e < graphLifted_.numberOfEdges(); ++e) {
		if(finalEdges[e]){
			size_t v0=graphLifted_.vertexOfEdge(e,0);
			size_t v1=graphLifted_.vertexOfEdge(e,1);
			if(liftedEdges[v0].count(v1)>0){
				tempGraphLifted.insertEdge(v0,v1);
				newLiftedCosts.push_back(liftedEdgeScore[e]);
			}
			else{
				auto edgeTest=graph_.findEdge(v0,v1);
				if(edgeTest.first){
					edgeScore[edgeTest.second]+=liftedEdgeScore[e];  //Compensate that the lifted edge has been removed
				}
				else{
					throw std::runtime_error("Error in lifted graph sparsification.");
				}

			}
		}
	}






	liftedEdgeScore=newLiftedCosts;

	graphLifted_=tempGraphLifted;
	std::cout<<"Left "<<newLiftedCosts.size()<<" lifted edges."<<std::endl;
	parameters.infoFile()<<"Left "<<newLiftedCosts.size()<<" lifted edges."<<std::endl;

	if(graphLifted_.numberOfEdges()!=newLiftedCosts.size()){
		std::cout<<"lifted edge number mismatch, lifted graph: "<<graphLifted_.numberOfEdges()<<", cost vector "<<newLiftedCosts.size()<<std::endl;
		parameters.infoFile()<<"lifted edge number mismatch, lifted graph: "<<graphLifted_.numberOfEdges()<<", cost vector "<<newLiftedCosts.size()<<std::endl;
	}
	else{
		std::cout<<"lifted edge number and lifted graph size match "<<std::endl;
		parameters.infoFile()<<"lifted edge number and lifted graph size match "<<std::endl;

	}
	parameters.infoFile().flush();

}





template<class T = size_t>
struct Data {

public:
	Data(DisjointStructure<> & problemGraph) :
		parameters(problemGraph.parameters) {

		//std::cout<<"data new constructor"<<std::endl;

		pGraph=problemGraph.getPGraph();
		pGraphLifted=problemGraph.getPGraphLifted();

		numberOfEdges = problemGraph.getGraph().numberOfEdges();
		numberOfVertices = problemGraph.getGraph().numberOfVertices();
		numberOfLiftedEdges = problemGraph.getGraphLifted().numberOfEdges();
		s=problemGraph.getSourceNode();
		t=problemGraph.getTerminalNode();

		costs=std::vector<double>();
		costs.insert(costs.begin(),problemGraph.getVerticesScore().begin(),problemGraph.getVerticesScore().end());
		costs.insert(costs.end(),problemGraph.getEdgesScore().begin(),problemGraph.getEdgesScore().end());
		costs.insert(costs.end(),problemGraph.getLiftedEdgesScore().begin(),problemGraph.getLiftedEdgesScore().end());


		//pReachable=problemGraph.getPReachable();
		pReachable=problemGraph.getPReachableNew();


		//if(parameters.isRestrictFrames()||parameters.isSparsify()){
			useTimeFrames=true;
			pVertexGroups=problemGraph.getPVertexGroups();
//	    }
//		else{
//			pVertexGroups=0;
//			useTimeFrames=false;
//		}

		pGraphComplete=0;
		pCompleteScore=0;

		useTimeTracks=false;

	}


	Data(CompleteStructure<>& cs,ConfigDisjoint<>& inpParameters) :
		parameters(inpParameters) {

		std::cout<<"new constructor of data"<<std::endl;
		parameters.infoFile()<<"new constructor of data"<<std::endl;
		parameters.infoFile().flush();
		pGraph=0;
		pGraphLifted=0;

		pGraphComplete=&cs.completeGraph;
		pCompleteScore=&cs.completeScore;

		pVertexGroups=&cs.vg;

		numberOfVertices=pGraphComplete->numberOfVertices();

		useTimeFrames=false;
		useTimeTracks=false;

		//numberOfVertices=pGraphComplete->numberOfVertices();
		std::cout<<"init pointers set"<<std::endl;
		parameters.infoFile()<<"init pointers set"<<std::endl;
		s=0;
		t=0;
		numberOfVertices=pGraphComplete->numberOfVertices();
		numberOfEdges=pGraphComplete->numberOfEdges();
		numberOfLiftedEdges=pGraphComplete->numberOfEdges();
		pReachable=0;
		//prepareGraphFromIntervalsDense(paths);


	}

	~Data(){
		std::cout<<"call data destructor "<<std::endl;
		parameters.infoFile()<<"call data destructor "<<std::endl;
		parameters.infoFile().flush();
		deleteCompleteGraph();
	}





	const andres::graph::Digraph<>& getGraph() {  //Consider removing and replacing with pointer getter
		return *pGraph;
	}

	 const andres::graph::Digraph<>& getGraphLifted()  { //Consider removing and replacing with pointer getter
		return *pGraphLifted;
	 }

	 const andres::graph::Digraph<>* getPGraph() {
		 return pGraph;
	 }

	 const andres::graph::Digraph<>* getPGraphLifted()  {
		 return pGraphLifted;
	 }



	 const ConfigDisjoint<>* getParameters()  {
	 		 return &parameters;
	 }


	size_t getSourceNode() const {
		return s;
	}


    size_t getTerminalNode()  {
		return t;
	}



    bool isReachable (size_t v,size_t w) {
    	if(v==s||w==t){
    		return true;
    	}
    	else if(w==s||v==t){
    		return false;
    	}
    	else if(pReachable==0){
    		return true;
    	}
    	else{
    		//return (*pReachable)[v][w];
    		return (*pReachable)[v].count(w)>0;
    	}
    }

	bool useTime () {
			return useTimeFrames;
	}

	bool useTimeForTracks () {
		return useTimeTracks;
	}


	size_t getGroupIndex(size_t v) {
		if(useTimeFrames){
			return pVertexGroups->getGroupIndex(v);
		}
		else throw std::runtime_error(std::string("Acces to time frames not enabled in this instance"));
	}


	VertexGroups<size_t>* getPVertexGroups(){
		if(useTimeFrames){
			return pVertexGroups;
		}
		else throw std::runtime_error(std::string("Acces to time frames not enabled in this instance"));
	}



	size_t getVertexVarIndex(size_t vertexIndex) const {
		assert(vertexIndex >= 0 && vertexIndex < numberOfVertices);
		return vertexIndex;
	}

	size_t getEdgeVarIndex(size_t edgeIndex) const {
		assert(edgeIndex >= 0 && edgeIndex < numberOfEdges);
		return edgeIndex + numberOfVertices;

	}

	size_t getLiftedEdgeVarIndex(size_t edgeIndex) const {
		assert(edgeIndex >= 0 && edgeIndex < numberOfLiftedEdges);
		return edgeIndex + numberOfVertices + numberOfEdges;
	}

	std::pair<std::string,size_t> getVariableType(size_t index) {
		assert(index >= 0 && index < costs.size());


		std::vector<size_t>breakPoints= {0,numberOfVertices,numberOfVertices+numberOfEdges,numberOfVertices+numberOfEdges+numberOfLiftedEdges};
		std::vector<std::string> type={"v-","e-","le-"};
		std::string outputType;
		size_t outIndex=0;
		for (int i = 0; i < breakPoints.size()-1; ++i) {
            if(index<breakPoints[i+1]){
            	outputType=type[i];
            	outIndex=index-breakPoints[i];
            	break;
            }
		}
		std::pair<std::string,size_t> output(outputType,outIndex);
		return output;

	}



	void graphFromIntervalsDense();
	void onePathGraph(std::vector<size_t>& path,bool zeroInOut=false);
	std::vector<size_t> cutsInOnePathSolution(std::vector<double>& labels);
	std::vector<std::vector<size_t>> pathsFromOnePathSolution(std::vector<double>& labels);
	bool smallPathsFromTimeBreaks(std::vector<std::vector<size_t>>& paths,bool finalCheckOnly=false);
	bool prepareGraphFromIntervalsDense(std::vector<std::vector<size_t>>& paths,bool finalCheckOnly=false);
	bool graphFromOptimizedPaths(std::vector<std::vector<size_t>>& paths,std::unordered_map<size_t,std::vector<size_t>>& breaks,bool finalCheckOnly);
	std::unordered_map<size_t,std::vector<size_t>> findTrajectoryBreaks(std::vector<std::vector<size_t>>& paths);
	//bool graphFromOptimizedPaths(std::vector<std::vector<size_t>>& paths,std::unordered_map<size_t,std::vector<size_t>>& breaks,bool finalCheckOnly);

	bool smallPathsFromTrajectoryBreaks(std::vector<std::vector<size_t>>& paths,bool finalCheckOnly);
	void prepareGraphFromPathsAndBreaks(std::vector<std::vector<size_t>>& paths,std::unordered_map<size_t,std::vector<size_t>>& breaks);
	std::vector<double> trackletsFromPaths(std::vector<std::vector<size_t>>& paths,std::unordered_map<size_t,std::set<size_t>>* timeBreaks=0);


	void outputSolution(std::vector<std::vector<size_t>>& paths,bool isIntervals=false);

	std::unordered_map<size_t,std::set<size_t>> findTimeBreaks(std::vector<std::vector<size_t>>& paths,bool origOnly=false);
	std::pair<double,double> evaluate(std::vector<std::vector<size_t>>& paths);
	std::vector<std::vector<size_t>> pathsFromSolution(std::vector<double>& intervalLabels,bool decodeTracklets=false,size_t shift=0);

	void readCompleteGraph();
	void deleteCompleteGraph();

	const std::vector<double> * getPCosts(){
		return &costs;
	}

	const std::vector<double>& getCosts(){
		return costs;
	}

	char getTask() const {
		return task;
	}

	const std::map<size_t, std::unordered_set<size_t> >& getTimeToTrGroup() const {
		if(useTimeTracks){
			return timeToTrGroup;
		}
		else{
			throw std::runtime_error("tracklet time groups now available in this settings");
		}
	}

	const std::vector<std::unordered_set<size_t> >& getTrackletsToGroups() const {
		if(useTimeTracks){
			return trackletsToGroups;
		}
		else{
			throw std::runtime_error("tracklet time groups now available in this settings");
		}
	}

	const std::vector<double>& getInitSolution() const {
		return initSolution;
	}

	size_t numberOfVertices;

	size_t numberOfEdges;

	size_t numberOfLiftedEdges;


	levinkov::Timer timer;

	disjointPaths::ConfigDisjoint<>& parameters;

private:
//	bool isReachableInTemp(size_t v,size_t w){
//		if(v==s||w==t){
//		    return true;
//		}
//		else if(w==s||v==t){
//			return false;
//		}
//		else if(pReachable==0){
//			return true;
//		}
//		else{
//			return (*pReachable)[v][w];
//		}
//
//	}
	void addSmallerPaths(std::vector<size_t>& newSmallPath);

	size_t s;
	size_t t;

	bool useTimeFrames;
	bool useTimeTracks;

	andres::graph::Digraph<> trackletGraph;
	andres::graph::Digraph<> trackletGraphLifted;

	andres::graph::Digraph<>* pGraph;
	andres::graph::Digraph<>* pGraphLifted;

	andres::graph::Digraph<>* pGraphComplete;

	VertexGroups<>* pVertexGroups;

	std::vector<double> costs;

	std::map<size_t,std::unordered_set<size_t>> timeToTrGroup;
	std::vector<std::unordered_set<size_t>> trackletsToGroups;

	//std::vector<std::vector<bool>>* pReachable;
	std::vector<std::unordered_set<size_t>>* pReachable;

	//std::vector<std::vector<bool>> newReachable;
	std::vector<std::unordered_set<size_t>> newReachable;
	std::vector<double>* pCompleteScore;

	std::vector<std::vector<size_t>> smallPaths;
	std::vector<int> vertexToPath;
	std::vector<std::vector<size_t>> pathsForInit;
	std::vector<double> initSolution;


	char task='C';

};




template<class T>
inline std::vector<double> Data<T>::trackletsFromPaths(std::vector<std::vector<size_t>>& paths,std::unordered_map<size_t,std::set<size_t>>* timeBreaks){
	std::cout<<"tracklets from paths"<<std::endl;
	parameters.infoFile()<<"tracklets from paths"<<std::endl;
	//parameters.output("tracklets from paths\n");

	vertexToPath=std::vector<int> (pGraphComplete->numberOfVertices(),-1);
	smallPaths=std::vector<std::vector<size_t>>();
	std::cout<<"find small paths"<<std::endl;
	parameters.infoFile()<<"find small paths"<<std::endl;
	parameters.infoFile().flush();
	//parameters.output("find small paths\n");
	size_t maxTime=0;
	std::vector<double> nodeCosts;

	bool createNodeCosts=parameters.getSmallIntervals()==0&&!parameters.isDenseTracklets();

	//bool halfIntervals=parameters.isOverlappingIntervals();
	bool halfIntervals=false;


	std::vector<size_t> vertexToInterval(pGraphComplete->numberOfVertices());
	//TODO s and t separately
	if(parameters.getSmallIntervals()>0){
		if(!halfIntervals){
			size_t I=parameters.getSmallIntervals();
			size_t trSize=parameters.getTrackletSize();
			size_t trNumber=I/trSize;
			if(I%trSize!=0){
				trNumber++;
			}
			std::vector<size_t> timeToTr(pVertexGroups->getMaxTime());

			for (int v = 0; v < pGraphComplete->numberOfVertices(); ++v) {
				size_t g=pVertexGroups->getGroupIndex(v);
				g--;  //compensate numbering from 1
				vertexToInterval[v]=trNumber*(g/I)+(g%I)/trSize+1;
				timeToTr[g]=vertexToInterval[v];
			}
//			std::cout<<"time to tr. number "<<std::endl;
//			for (int t = 0; t <= pVertexGroups->getMaxTime(); ++t) {
//			    std::cout<<t<<":"<<timeToTr[t]<<std::endl;
//			}
			std::cout<<"vertex to interval done"<<std::endl;
			parameters.infoFile()<<"vertex to interval done"<<std::endl;
			parameters.infoFile().flush();
		}
//		else{
//			size_t I=parameters.getSmallIntervals()/2;
//			size_t trSize=parameters.getTrackletSize();
//			size_t trNumber=I/trSize;
//			if(I%trSize!=0){
//				trNumber++;
//			}
//			for (int v = 0; v < numberOfVertices-2; ++v) {
//				size_t g=pVertexGroups->getGroupIndex(v);
//				g--;  //compensate numbering from 1
//				size_t indexI=0;
//				size_t indexTr=0;
//				size_t finalIndex=0;
//				if(g<I+I/2){
//					indexI=0;
//					indexTr=g/trSize;
//					finalIndex=indexTr+1;
//				}
//				else{
//					g=g-I/2;
//					indexI=g/I;
//					indexTr=(g-indexI*I)/trSize;
//					finalIndex=(I+I/2)/trSize+2+indexI*(trNumber)+indexTr;
//				}
//				vertexToInterval[v]=finalIndex;
//			}
//		}
	}
	else if(parameters.getTrackletSize()>0){
		size_t trSize=parameters.getTrackletSize();
		for (int v = 0; v < numberOfVertices; ++v) {
			size_t g=pVertexGroups->getGroupIndex(v);
			g--;  //compensate numbering from 1
			vertexToInterval[v]=g/trSize+1;
		}
	}

	for (int i = 0; i < paths.size(); ++i) {
	//	std::cout<<"path "<<i<<std::endl;
		double currentNodeCost=0;
		std::vector<size_t> path;
		for (int j = 0; j < paths[i].size(); ++j) {
			size_t v=paths[i][j];
			path.push_back(v);
			//vertexToPath[v]=smallPaths.size();

			if(j!=paths[i].size()-1){
				size_t w=paths[i][j+1];


				if(vertexToInterval[v]!=vertexToInterval[w]||(*timeBreaks)[i].count(pVertexGroups->getGroupIndex(v))>0){
					smallPaths.push_back(path);
					nodeCosts.push_back(currentNodeCost);
					path=std::vector<size_t>();
					currentNodeCost=0;

				}
				else if(createNodeCosts){
					auto findEdge=pGraph->findEdge(v,w);
					currentNodeCost+=costs[getEdgeVarIndex(findEdge.second)];
				}
			}
			else{
				smallPaths.push_back(path);
				nodeCosts.push_back(currentNodeCost);
				path=std::vector<size_t>();
									currentNodeCost=0;
			}


		}

	}

	std::cout<<"tracklets from paths"<<std::endl;
	parameters.infoFile()<<"tracklets from paths"<<std::endl;
	parameters.infoFile().flush();
	//parameters.output("small paths created. Creating single node paths\n");
	for (int i = 0; i < smallPaths.size(); ++i) {
		for (int j = 0; j < smallPaths[i].size(); ++j) {
			if(vertexToPath[smallPaths[i][j]]!=-1){
				std::cout<<"multiple occurences of vertex "<<smallPaths[i][j]<<std::endl;
				throw std::runtime_error("error in optimized paths, multiple occurrences of one vertex");
			}
			else{
				vertexToPath[smallPaths[i][j]]=i;
			}
		}
	}


	for (int i = 0; i < vertexToPath.size(); ++i) {
		if(vertexToPath[i]==-1){
			std::vector<size_t> path;
			path.push_back(i);
			vertexToPath[i]=smallPaths.size();
			smallPaths.push_back(path);
			nodeCosts.push_back(0);
		}
	}

    return nodeCosts;
}



template<class T>
inline void Data<T>::deleteCompleteGraph(){
	if(pGraphComplete!=0&&parameters.getSmallIntervals()==0){
		std::cout<<"delete complete graph "<<std::endl;
		parameters.infoFile()<<"delete complete graph "<<std::endl;
		parameters.infoFile().flush();
		//parameters.output("delete complete graph\n");
		delete pGraphComplete;
		delete pCompleteScore;
	}
}

//TODO make one procedure and use it in CS too
template<class T>
inline void Data<T>::readCompleteGraph(){
	if(pGraphComplete==0){

		std::string line;
		std::ifstream data;
		try{
			data.open(parameters.getGraphFileName());
			if(!data){

				throw std::system_error(errno, std::system_category(), "failed to open "+parameters.getGraphFileName());

			}
			char delim=',';
			std::getline(data, line);

			//size_t origVertexNumber=std::stoul(line);
			size_t origVertexNumber=pVertexGroups->getMaxVertex()+1;
			pGraphComplete=new andres::graph::Digraph<>(origVertexNumber);
			pCompleteScore=new std::vector<double>();
			double objValue=0;

			std::cout << "Read complete graph" << std::endl;
			parameters.infoFile()<< "Read complete graph" << std::endl;
			std::vector<std::string> strings;


			std::cout<<"Skipping vertices in file. "<<std::endl;
			parameters.infoFile()<<"Skipping vertices in file. "<<std::endl;
			parameters.infoFile().flush();
			//Vertices that are not found have score=0. Appearance and disappearance cost are read here.
			while (std::getline(data, line) && !line.empty()) {

			}

			std::cout<<"Reading base edges from file. "<<std::endl;
			parameters.infoFile()<<"Reading base edges from file. "<<std::endl;
			parameters.infoFile().flush();
			size_t maxLabel=0;
			while (std::getline(data, line) && !line.empty()) {

				strings = split(line, delim);

				unsigned int v = std::stoul(strings[0]);
				unsigned int w = std::stoul(strings[1]);
				size_t l0=pVertexGroups->getGroupIndex(v);
				size_t l1=pVertexGroups->getGroupIndex(w);
				if(v>=origVertexNumber||w>=origVertexNumber) continue;
				if(l1-l0<=parameters.getMaxTimeGapComplete()){
					double score = std::stod(strings[2]);
					pGraphComplete->insertEdge(v,w);
					pCompleteScore->push_back(score);
				}

			}

			data.close();
		}
		catch (std::system_error& er) {
			std::clog << er.what() << " (" << er.code() << ")" << std::endl;

		}
	}

}




template<class T>
inline void Data<T>::onePathGraph(std::vector<size_t>& path,bool zeroInOut){
	//std::cout<<"call one path graph"<<std::endl;
	readCompleteGraph();
	std::vector<size_t> newToOrigVertices(path.size());
	//andres::graph::Digraph<> oneTrackGraph(path.size()+2);
	//andres::graph::Digraph<> oneTrackGraphLifted(path.size()+2);
	pGraph=new andres::graph::Digraph<>(path.size()+2);
	pGraphLifted=new andres::graph::Digraph<>(path.size()+2);
	std::vector<double> newBaseCosts;
	std::vector<double> newLiftedCosts;

	newToOrigVertices=path;
	costs=std::vector<double>(path.size()+2);
	numberOfVertices=path.size()+2;
	newReachable=std::vector<std::unordered_set<size_t>>(path.size()+2);
	s=path.size();
	t=path.size()+1;

	//std::cout<<"construct edges"<<std::endl;


	for (int i = 0; i < path.size()-1; ++i) {
		size_t v=path[i];
		newReachable[i].insert(i);
		newReachable[s].insert(i);
		newReachable[i].insert(t);
		for (int j = i+1; j < path.size(); ++j) {
			size_t w=path[j];
			newReachable[i].insert(j);
			auto findEdge=pGraphComplete->findEdge(v,w);
			if(findEdge.first){
				pGraph->insertEdge(i,j);
				pGraphLifted->insertEdge(i,j);
				double cost=(*pCompleteScore)[findEdge.second];
				newBaseCosts.push_back(cost);
				size_t t0=pVertexGroups->getGroupIndex(v);
				size_t t1=pVertexGroups->getGroupIndex(w);
				if(cost<=0||t1-t0<=parameters.getRepulsiveTimeGap()){
					if(t1-t0<=parameters.getTrTimeGapLifted()){
						newLiftedCosts.push_back(cost);
					}
					else{
						newLiftedCosts.push_back(0);
					}
				}
				else{
					newLiftedCosts.push_back(0);
				}

			}

		}
	}


	//std::cout<<"insert s t edges"<<std::endl;
	for (int i = 0; i < path.size(); ++i) {
		//std::cout<<"vertex "<<v;
		pGraph->insertEdge(s,i);
		pGraph->insertEdge(i,t);
		if(zeroInOut){
			newBaseCosts.push_back(0);
			newBaseCosts.push_back(0);
			//std::cout<<"zero cost inserted "<<std::endl;
		}
		else{
			newBaseCosts.push_back(parameters.getInputCost());
			newBaseCosts.push_back(parameters.getOutputCost());
			//std::cout<<"non zero cost inserted "<<std::endl;
		}

	}


	//std::cout<<"set pointers"<<std::endl;
//	pGraph=new andres::graph::Digraph<>(oneTrackGraph);
//	pGraphLifted=new andres::graph::Digraph<>(oneTrackGraphLifted);
	costs.insert(costs.end(),newBaseCosts.begin(),newBaseCosts.end());
	costs.insert(costs.end(),newLiftedCosts.begin(),newLiftedCosts.end());
	pReachable=&newReachable;
	useTimeFrames=false;
	useTimeTracks=false;
	numberOfEdges=pGraph->numberOfEdges();
	numberOfLiftedEdges=pGraphLifted->numberOfEdges();

}




template<class T>
inline std::vector<size_t> Data<T>::cutsInOnePathSolution(std::vector<double>& labels){
	std::vector<size_t> cutIndices;
	for (size_t v = 0; v < numberOfVertices-3; ++v) {
		size_t w=v+1;
		auto findEdge=pGraph->findEdge(v,w);
		if(!findEdge.first){
			throw std::runtime_error("error in one path graph solution");
		}
		else{
			if(labels[getEdgeVarIndex(findEdge.second)]<0.5){
				cutIndices.push_back(v);
			}

		}
	}

	delete pGraph;
	delete pGraphLifted;
	pGraph=0;
	pGraphLifted=0;

	return cutIndices;
}



template<class T>
inline bool Data<T>::smallPathsFromTrajectoryBreaks(std::vector<std::vector<size_t>>& paths,bool finalCheckOnly){
	std::unordered_map<size_t,std::vector<size_t>> trajBreaks=findTrajectoryBreaks(paths);
	numberOfVertices=pGraphComplete->numberOfVertices();

	if(!finalCheckOnly||!trajBreaks.empty()){

		prepareGraphFromPathsAndBreaks(paths,trajBreaks);


		std::cout<<"evaluate small paths"<<std::endl;
		parameters.infoFile()<<"evaluate small paths"<<std::endl;
		parameters.infoFile().flush();
		trajBreaks=findTrajectoryBreaks(smallPaths);

		while(!trajBreaks.empty()){
			std::cout<<"cutting small paths "<<std::endl;
			parameters.infoFile()<<"cutting small paths "<<std::endl;
			parameters.infoFile().flush();
			std::vector<std::vector<size_t>> oldSmallPaths=smallPaths;
			numberOfVertices=pGraphComplete->numberOfVertices();
			prepareGraphFromPathsAndBreaks(oldSmallPaths,trajBreaks);
			trajBreaks=findTrajectoryBreaks(smallPaths);
		}

		return true;
	}
	else{
		return false;
	}


}



template<class T>
inline bool Data<T>::smallPathsFromTimeBreaks(std::vector<std::vector<size_t>>& paths,bool finalCheckOnly){
	std::unordered_map<size_t,std::set<size_t>> timeBreaks=findTimeBreaks(paths,true);
	numberOfVertices=pGraphComplete->numberOfVertices();

	if(!finalCheckOnly||!timeBreaks.empty()){

		trackletsFromPaths(paths,&timeBreaks);

		std::cout<<"evaluate small paths"<<std::endl;
		parameters.infoFile()<<"evaluate small paths"<<std::endl;
		parameters.infoFile().flush();
		timeBreaks=findTimeBreaks(smallPaths,true);

		while(!timeBreaks.empty()){
			std::cout<<"cutting small paths "<<std::endl;
			parameters.infoFile()<<"cutting small paths "<<std::endl;
			parameters.infoFile().flush();
			std::vector<std::vector<size_t>> oldSmallPaths=smallPaths;
			numberOfVertices=pGraphComplete->numberOfVertices();
			trackletsFromPaths(oldSmallPaths,&timeBreaks);
			timeBreaks=findTimeBreaks(smallPaths,true);
		}

		return true;
	}
	else{
		return false;
	}


}


template<class T>
inline void Data<T>::addSmallerPaths(std::vector<size_t>& newSmallPath){
	if(parameters.getTrackletSize()>0){
		std::vector<size_t> initPath;
		size_t vFirst=newSmallPath[0];
		size_t vLast=(*newSmallPath.rbegin());
		size_t tFirst=pVertexGroups->getGroupIndex(vFirst);
		size_t tLast=pVertexGroups->getGroupIndex(vLast);
		if((tFirst-1)/parameters.getTrackletSize()!=(tLast-1)/parameters.getTrackletSize()){
			size_t currentInterval=(tFirst-1)/parameters.getTrackletSize();
			std::vector<size_t> newSmallerPath;
			for(size_t v:newSmallPath){
				size_t newInterval=(pVertexGroups->getGroupIndex(v)-1)/parameters.getTrackletSize();
				if(newInterval==currentInterval){
					newSmallerPath.push_back(v);
				}
				else{
					initPath.push_back(smallPaths.size());
					smallPaths.push_back(newSmallerPath);
					newSmallerPath=std::vector<size_t>();
					newSmallerPath.push_back(v);
					currentInterval=newInterval;
				}

			}
			initPath.push_back(smallPaths.size());
			smallPaths.push_back(newSmallerPath);

		}
		else{
			initPath.push_back(smallPaths.size());
			smallPaths.push_back(newSmallPath);
		}
		pathsForInit.push_back(initPath);
	}
	else{
		smallPaths.push_back(newSmallPath);
	}

}




template<class T>
inline void Data<T>::prepareGraphFromPathsAndBreaks(std::vector<std::vector<size_t>>& paths,std::unordered_map<size_t,std::vector<size_t>>& breaks){


	std::cout<<"tracklets from optimized paths"<<std::endl;
	parameters.infoFile()<<"tracklets from optimized paths"<<std::endl;
	//parameters.output("tracklets from paths\n");

	vertexToPath=std::vector<int> (pGraphComplete->numberOfVertices(),-1);
	smallPaths=std::vector<std::vector<size_t>>();
	pathsForInit=std::vector<std::vector<size_t>>();

//	std::cout<<"all paths"<<std::endl;
//	for(auto& path:paths){
//		for(size_t v:path){
//			std::cout<<v<<", ";
//		}
//		std::cout<<std::endl;
//	}
//	std::cout<<"display breaks"<<std::endl;
//	for(auto& pair:breaks){
//		std::cout<<pair.first<<": ";
//		for(size_t br:pair.second){
//			std::cout<<br<<", ";
//		}
//		std::cout<<std::endl;
//	}


	int counter=1;
	for (int i = 0; i < paths.size(); ++i) {
		if(breaks.count(i)>0&&!breaks[i].empty()){
			//updateDone=true;
			size_t startTrack=0;
			for (int j = 0; j < breaks[i].size(); ++j) {
				if(breaks[i][j]>startTrack){
					std::vector<size_t> newSmallPath(paths[i].begin()+startTrack,paths[i].begin()+breaks[i][j]+1);
					addSmallerPaths(newSmallPath);
				}
				startTrack=breaks[i][j]+1;
			}
			if(startTrack<paths[i].size()-1){
				std::vector<size_t> newSmallPath(paths[i].begin()+startTrack,paths[i].end());
				addSmallerPaths(newSmallPath);
			}

		}
		else{
			addSmallerPaths(paths[i]);
		}
	}


	std::cout<<"display small paths"<<std::endl;
	for (int i = 0; i < smallPaths.size(); ++i) {
		for (int j = 0; j < smallPaths[i].size(); ++j) {
			std::cout<<smallPaths[i][j]<<", ";
			if(vertexToPath[smallPaths[i][j]]!=-1){
				std::cout<<"multiple occurences of vertex "<<smallPaths[i][j]<<std::endl;
				throw std::runtime_error("error in optimized paths, multiple occurrences of one vertex");
			}
			else{
				vertexToPath[smallPaths[i][j]]=i;
			}
		}
		std::cout<<std::endl;
	}
	for (int i = 0; i < vertexToPath.size(); ++i) {
		if(vertexToPath[i]==-1){
			std::vector<size_t> oneVertexPath;
			oneVertexPath.push_back(i);
			vertexToPath[i]=smallPaths.size();
			smallPaths.push_back(oneVertexPath);
		}
	}

}





template<class T>
inline bool Data<T>::graphFromOptimizedPaths(std::vector<std::vector<size_t>>& paths,std::unordered_map<size_t,std::vector<size_t>>& breaks,bool finalCheckOnly){


	if(finalCheckOnly&&breaks.empty()){
		return false;
	}
	else{
		prepareGraphFromPathsAndBreaks(paths,breaks);


		graphFromIntervalsDense();
		//return updateDone;
		return true;

	}
}


template<class T>
inline bool Data<T>::prepareGraphFromIntervalsDense(std::vector<std::vector<size_t>>& paths,bool finalCheckOnly){

	readCompleteGraph();

	//bool createGraph=smallPathsFromTimeBreaks(paths,finalCheckOnly);
	bool createGraph=smallPathsFromTrajectoryBreaks(paths,finalCheckOnly);

	if(createGraph){
		graphFromIntervalsDense();
		return true;
	}
	else
	{
		return false;

	}
}


template<class T>
inline void Data<T>::graphFromIntervalsDense(){

	    initSolution=std::vector<double>();
	    std::cout<<"all small paths"<<std::endl;
	    std::ofstream pathsFile;
	    pathsFile.open(parameters.getOutputFileName()+"paths.txt",std::ofstream::out | std::ofstream::app);

	    for (int i = 0; i < vertexToPath.size(); ++i) {
			size_t p=vertexToPath[i];
			if(smallPaths[p][0]==i){
				for (int j = 0; j < smallPaths[p].size(); ++j) {
					pathsFile<<smallPaths[p][j]<<",";
				}
				pathsFile<<std::endl;
			}
		}
	    pathsFile<<std::endl<<std::endl;

	    pathsFile.close();


		andres::graph::Digraph<> trGraphBase(smallPaths.size()+2);//Maybe without +2
		andres::graph::Digraph<> trGraphLifted(smallPaths.size()+2);
		std::vector<double> nodeCosts(smallPaths.size()+2,0);

		std::cout<<"small paths size "<<smallPaths.size()<<std::endl;
		parameters.infoFile()<<"small paths size "<<smallPaths.size()<<std::endl;

		std::unordered_map<size_t,size_t> insideBaseEdges; //would be enough to have only simple map
		std::unordered_map<size_t,std::unordered_set<size_t>> originalEdges; //would be enough to have map to set


		std::cout<<"check time consistency of tracklets"<<std::endl;
		parameters.infoFile()<<"check time consistency of tracklets"<<std::endl;
		parameters.infoFile().flush();


		for (int i = 0; i < smallPaths.size(); ++i) {
			size_t firstVertex=*(smallPaths[i].begin());
			size_t lastVertex=*(smallPaths[i].rbegin());
			size_t gFirst=vertexToPath[firstVertex];
			size_t gLast=vertexToPath[lastVertex];
			if(gFirst!=i||gLast!=i){
				std::cout<<"gFirst "<<gFirst<<", gLast "<<gLast<<", i "<<i<<std::endl;
				parameters.infoFile()<<"Error: Mismatch in tracklet numbers and vertexTOPath array"<<std::endl;
				parameters.infoFile()<<"gFirst "<<gFirst<<", gLast "<<gLast<<", i "<<i<<std::endl;
				parameters.infoFile().flush();
				throw std::runtime_error("Mismatch in tracklet numbers and vertexTOPath array");
			}

			size_t firstTime=pVertexGroups->getGroupIndex(firstVertex);
			size_t lastTime=pVertexGroups->getGroupIndex(lastVertex);

			//find new edges
			for (int j = 0; j < i; ++j) {
				size_t firstVertex2=*(smallPaths[j].begin());
				size_t lastVertex2=*(smallPaths[j].rbegin());
				size_t firstTime2=pVertexGroups->getGroupIndex(firstVertex2);
				size_t lastTime2=pVertexGroups->getGroupIndex(lastVertex2);

				size_t gFirst2=vertexToPath[firstVertex2];
				size_t gLast2=vertexToPath[lastVertex2];
				if(gFirst2!=j||gLast2!=j){
					parameters.infoFile()<<"Error: Mismatch in tracklet numbers and vertexTOPath array"<<std::endl;
					parameters.infoFile()<<"gFirst "<<gFirst<<", gLast "<<gLast<<", i "<<i<<std::endl;
					parameters.infoFile().flush();
					std::cout<<"gFirst2 "<<gFirst2<<", gLast2 "<<gLast2<<", j "<<j<<std::endl;
					throw std::runtime_error("Mismatch in tracklet numbers and vertexTOPath array");
				}


				if(lastTime<firstTime2){
					if(firstTime2-lastTime<=parameters.getTrTimeGapBase()){
						originalEdges[lastVertex].insert(firstVertex2);
						trGraphBase.insertEdge(i,j);
					}
					if(firstTime2-lastTime<=parameters.getTrTimeGapLifted()){
						trGraphLifted.insertEdge(i,j);
					}
				}
				else if(lastTime2<firstTime){
					if(firstTime-lastTime2<=parameters.getTrTimeGapBase()){
						originalEdges[lastVertex2].insert(firstVertex);
						trGraphBase.insertEdge(j,i);
					}
					if(firstTime-lastTime2<=parameters.getTrTimeGapLifted()){
						trGraphLifted.insertEdge(j,i);
					}

				}
			}
			//Will be used for node costs
			for (int j = 0; j < smallPaths[i].size()-1; ++j) {
				size_t v=smallPaths[i][j];
				size_t w=smallPaths[i][j+1];
				insideBaseEdges[v]=w;
			}
		}

		std::cout<<"base for graphs, edges: "<<trGraphBase.numberOfEdges()<<", lifted edges: "<<trGraphLifted.numberOfEdges()<<std::endl;
		parameters.infoFile()<<"base for graphs, edges: "<<trGraphBase.numberOfEdges()<<", lifted edges: "<<trGraphLifted.numberOfEdges()<<std::endl;



		std::vector<double> liftedEdgeCosts(trGraphLifted.numberOfEdges());
		std::vector<double> baseEdgeCosts(trGraphBase.numberOfEdges());


		std::unordered_map<size_t,std::unordered_set<size_t>> usedEdges;
		std::cout<<"insert original costs into tracklet costs"<<std::endl;
		parameters.infoFile()<<"insert original costs into tracklet costs"<<std::endl;
		for (int e = 0; e < pGraphComplete->numberOfEdges(); ++e) {

			size_t v=pGraphComplete->vertexOfEdge(e,0);
			size_t w=pGraphComplete->vertexOfEdge(e,1);

			int g0=vertexToPath[v];
			int g1=vertexToPath[w];

			usedEdges[g0].insert(g1);

			double cost=(*pCompleteScore)[e];
			double liftedCost=cost;
			size_t t0=pVertexGroups->getGroupIndex(v);
			size_t t1=pVertexGroups->getGroupIndex(w);
			if(t1-t0>parameters.getRepulsiveTimeGap()&&cost>0){
				liftedCost=0;
			}
			if(t1-t0>parameters.getTrTimeGapLifted()){
				liftedCost=0;
			}



			if(g0!=-1&&g1!=-1){
				if(g0==g1){
					if(insideBaseEdges.count(v)>0&&insideBaseEdges[v]==w){
						nodeCosts[g0]+=cost+liftedCost;
					}
					else{
						nodeCosts[g0]+=liftedCost;
					}

				}
				else{
					auto findEdge1=trGraphBase.findEdge(g0,g1);
					auto findEdge2=trGraphLifted.findEdge(g0,g1);

					if(findEdge1.first){
						if(originalEdges.count(v)>0&&originalEdges[v].count(w)>0){
//							if(originalEdges[v][w]!=findEdge1.second){
//								std::cout<<"orig edges "<<originalEdges[v][w]<<", edge index "<<findEdge1.second<<std::endl;
//								throw std::runtime_error("Error in building tracklet graph base");
//							}
							baseEdgeCosts[findEdge1.second]=cost;
						}


					}
					if(findEdge2.first){
					//	if(cost<=parameters.getNegativeThresholdLifted()||cost>=parameters.getPositiveThresholdLifted()){
							liftedEdgeCosts[findEdge2.second]+=liftedCost;
					//	}
					}

				}

			}
			else{
				throw std::runtime_error("Unassigned vertex for tracklet graph.");
			}

		}


		numberOfVertices=smallPaths.size()+2;
		trackletGraph=andres::graph::Digraph<>(smallPaths.size()+2);
		trackletGraphLifted=andres::graph::Digraph<>(smallPaths.size()+2);
		std::vector<double> baseEdgeCostNew;
		std::vector<double> liftedEdgeCostNew;
		double baseThreshold=parameters.getInputCost()+parameters.getOutputCost();

		std::cout<<"Construct base tracklet graph"<<std::endl;
		parameters.infoFile()<<"Construct base tracklet graph"<<std::endl;
		parameters.infoFile().flush();

		for (int e = 0; e < trGraphBase.numberOfEdges(); ++e) {
			size_t v0=trGraphBase.vertexOfEdge(e,0);
			size_t v1=trGraphBase.vertexOfEdge(e,1);
			auto findEdge=trGraphLifted.findEdge(v0,v1);
			double cost;
			if(findEdge.first){
				cost=baseEdgeCosts[e]+liftedEdgeCosts[findEdge.second];
			}
			else{
				cost=baseEdgeCosts[e];
			}

			if((parameters.isAllBaseTracklet()||cost<=baseThreshold)&&usedEdges[v0].count(v1)>0){
			//if(usedEdges[v0].count(v1)>0){
				trackletGraph.insertEdge(v0,v1);
				baseEdgeCostNew.push_back(baseEdgeCosts[e]);
				//baseEdgeCostNew.push_back(0);
			}
//			trackletGraph.insertEdge(v0,v1);
//			baseEdgeCostNew.push_back(0);



		}


		size_t newS=smallPaths.size();
		size_t newT=smallPaths.size()+1;
		s=newS;
		t=newT;

		std::cout<<"s t edges insert"<<std::endl;
		parameters.infoFile()<<"s t edges insert"<<std::endl;
		parameters.infoFile().flush();
		for (int v = 0; v < trackletGraph.numberOfVertices()-2; ++v) {
			trackletGraph.insertEdge(newS,v);
			baseEdgeCostNew.push_back(parameters.getInputCost());
			trackletGraph.insertEdge(v,newT);
			baseEdgeCostNew.push_back(parameters.getOutputCost());
		}


		//newReachable=initReachable(trackletGraph,parameters);
		newReachable=initReachableSet(trackletGraph,parameters);
		pReachable=&newReachable;



		std::cout<<"Find unnecessary lifted edges"<<std::endl;
		parameters.infoFile()<<"Find unnecessary lifted edges"<<std::endl;
		parameters.infoFile().flush();

		std::unordered_map<size_t,std::set<size_t>> leToKeep;
		for (int v = 0; v < trackletGraph.numberOfVertices()-2; ++v) {
			std::unordered_set<size_t> alternativePath;
			for (int i = 0; i < trackletGraph.numberOfEdgesFromVertex(v); ++i) {
				size_t w=trackletGraph.vertexFromVertex(v,i);
				for(size_t u:newReachable[w]){
					if(u!=w) alternativePath.insert(u);
				}
							}
			for (int i = 0; i < trGraphLifted.numberOfEdgesFromVertex(v); ++i) {
				size_t w=trGraphLifted.vertexFromVertex(v,i);
				if(w!=t){
					if(alternativePath.count(w)>0) leToKeep[v].insert(w);
				}

			}
		}


		std::cout<<"Construct tracklet graph lifted"<<std::endl;
		parameters.infoFile()<<"Construct tracklet graph lifted"<<std::endl;
		parameters.infoFile().flush();
		for (int e = 0; e < trGraphLifted.numberOfEdges(); ++e) {
			size_t v0=trGraphLifted.vertexOfEdge(e,0);
			size_t v1=trGraphLifted.vertexOfEdge(e,1);
			//if(usedEdges[v0].count(v1)>0&&newReachable[v0][v1]){
			if(usedEdges[v0].count(v1)>0&&newReachable[v0].count(v1)>0){
				if(leToKeep[v0].count(v1)>0){
					trackletGraphLifted.insertEdge(v0,v1);
					liftedEdgeCostNew.push_back(liftedEdgeCosts[e]);
				}
				else{
					auto findEdge=trackletGraph.findEdge(v0,v1);
					if(findEdge.first){
						baseEdgeCostNew[findEdge.second]+=liftedEdgeCosts[e];
					}
					else{
						throw std::runtime_error("Error in tracklet graph construction");
					}
				}

			}
		}





		//TODO order the paths according the intervals they belong to. One interval corresponds to a new time layer.



		//maxTime++;






		std::cout<<"set pointers"<<std::endl;
		parameters.infoFile()<<"set pointers"<<std::endl;
		parameters.infoFile().flush();
		pGraph=&trackletGraph;
		pGraphLifted=&trackletGraphLifted;
		numberOfEdges=pGraph->numberOfEdges();
		numberOfLiftedEdges=pGraphLifted->numberOfEdges();

		//costs=std::vector<double>(numberOfVertices,0);
		costs=nodeCosts;
		costs.insert(costs.end(),baseEdgeCostNew.begin(),baseEdgeCostNew.end());
		costs.insert(costs.end(),liftedEdgeCostNew.begin(),liftedEdgeCostNew.end());
		std::cout<<"number of vertices "<<numberOfVertices<<std::endl;
		std::cout<<"number of edges "<<numberOfEdges<<std::endl;
		std::cout<<"number of lifted edges "<<numberOfLiftedEdges<<std::endl;
		std::cout<<"cost size "<<costs.size()<<" sum "<<numberOfEdges+numberOfLiftedEdges+numberOfVertices<<std::endl;


		std::cout<<"paths for init "<<std::endl;
		if(task=='T'&&!pathsForInit.empty()){ //Do not do this in the first iteration where paths are obtained from intervals
			initSolution=std::vector<double>(costs.size(),0);
			for (int i = 0; i < pathsForInit.size(); ++i) {
				bool terminatedEarlier=false;
				//std::cout<<"path "<<i<<": ";
				auto findEdge=pGraph->findEdge(s,pathsForInit[i][0]);
				if(!findEdge.first){
					throw std::runtime_error("Edge from start vertex missing in tracklet graph");
				}
				else{
					initSolution[getEdgeVarIndex(findEdge.second)]=1;
				}
				for (int j = 0; j < pathsForInit[i].size()-1; ++j) {
					size_t v=pathsForInit[i][j];
					//std::cout<<v<<",";
					initSolution[getVertexVarIndex(v)]=1;

					size_t w=pathsForInit[i][j+1];
					findEdge=pGraph->findEdge(v,w);
					if(!findEdge.first){
						std::cout<<"Warning: Base graph does not contain possibly good edge "<<v<<","<<w<<std::endl;
						findEdge=pGraph->findEdge(v,t);
						if(!findEdge.first){
							throw std::runtime_error("Edge to end vertex missing in tracklet graph");
						}
						else{
							initSolution[getEdgeVarIndex(findEdge.second)]=1;
							terminatedEarlier=true;
						}
						break;
					}
					else{
						initSolution[getEdgeVarIndex(findEdge.second)]=1;
					}
					for (int k = 0; k <= j; ++k) {
						size_t u=pathsForInit[i][k];
						findEdge=pGraphLifted->findEdge(u,w);
						if(findEdge.first){
							initSolution[getLiftedEdgeVarIndex(findEdge.second)]=1;
						}
					}

				}
				if(!terminatedEarlier){
					size_t w=*(pathsForInit[i].rbegin());
					//std::cout<<w<<std::endl;
					initSolution[getVertexVarIndex(w)]=1;
					auto findEdge=pGraph->findEdge(w,t);
					if(!findEdge.first){
						throw std::runtime_error("Edge to end vertex missing in tracklet graph");
					}
					else{
						initSolution[getEdgeVarIndex(findEdge.second)]=1;
					}
				}


				//TODO last vertex to one and edge to t to one
			}



		}

		parameters.infoFile()<<"number of vertices "<<numberOfVertices<<std::endl;
		parameters.infoFile()<<"number of edges "<<numberOfEdges<<std::endl;
		parameters.infoFile()<<"number of lifted edges "<<numberOfLiftedEdges<<std::endl;
		parameters.infoFile()<<"cost size "<<costs.size()<<" sum "<<numberOfEdges+numberOfLiftedEdges+numberOfVertices<<std::endl;
		parameters.infoFile().flush();



		std::cout<<"new graphs complete"<<std::endl;
		parameters.infoFile()<<"new graphs complete"<<std::endl;
		parameters.infoFile().flush();
		useTimeFrames=false;




		trackletsToGroups=std::vector<std::unordered_set<size_t>> (smallPaths.size());
		timeToTrGroup.clear();
		size_t timeCounter=1;
		size_t maximalTime=pVertexGroups->getMaxTime();
		size_t pLast=1;
		const std::vector<size_t>& firstGroup= pVertexGroups->getGroupVertices(1);
		for (size_t v:firstGroup) {
			timeToTrGroup[1].insert(vertexToPath[v]);
		}


		timeCounter++;
		while(timeCounter<=maximalTime){
			const std::vector<size_t>& nextGroup= pVertexGroups->getGroupVertices(timeCounter);
			std::unordered_set<size_t> newGroup;
			for (size_t v:nextGroup) {
				newGroup.insert(vertexToPath[v]);
			}


			bool containedNewInOld=false;
			bool containedOldIneNew=false;
			if(newGroup.size()<=timeToTrGroup[pLast].size()){
				containedNewInOld=true;
				for(size_t tr:newGroup){
					if(timeToTrGroup[pLast].count(tr)==0){
						containedNewInOld=false;
						break;
					}
				}
				if(!containedNewInOld){
					for(size_t v:timeToTrGroup[pLast]){
						trackletsToGroups[v].insert(pLast);
					}
					timeToTrGroup[timeCounter]=newGroup;
					pLast=timeCounter;
				}
				//if(containedNewInOld) do nothing, else create new group

			}
			else if(newGroup.size()>timeToTrGroup[pLast].size()){
				containedOldIneNew=true;
				for(size_t tr:timeToTrGroup[pLast]){
					if(newGroup.count(tr)==0){
						containedOldIneNew=false;
						break;
					}
				}
				if(containedOldIneNew){
					timeToTrGroup[pLast]=newGroup;
				}
				else{
					for(size_t v:timeToTrGroup[pLast]){
						trackletsToGroups[v].insert(pLast);
					}
					timeToTrGroup[timeCounter]=newGroup;
					pLast=timeCounter;
				}
			}
			timeCounter++;
		}
		for(size_t v:timeToTrGroup[pLast]){
			trackletsToGroups[v].insert(pLast);
		}

		useTimeTracks=true;
		useTimeFrames=false;
		task='T';

}







template<class T>
inline std::vector<std::vector<size_t>>Data<T>::pathsFromSolution(std::vector<double>& intervalLabels, bool decodeTracklets,size_t shift){
	std::vector<std::vector<size_t>> paths;
	for (int i = 0; i < pGraph->numberOfEdgesFromVertex(s); ++i) {
		size_t edgeIndex=pGraph->edgeFromVertex(s,i);
		if(intervalLabels[getEdgeVarIndex(edgeIndex)]>0.5){
			std::vector<size_t> path;
			size_t v=pGraph->vertexFromVertex(s,i);
			while(v!=t){
				path.push_back(v+shift);
				for (int j = 0; j < pGraph->numberOfEdgesFromVertex(v); ++j) {
					size_t edge=pGraph->edgeFromVertex(v,j);
					if(intervalLabels[getEdgeVarIndex(edge)]>0.5){
						v=pGraph->vertexFromVertex(v,j);
						break;
					}
				}
			}
			paths.push_back(path);
		}
	}

	if(!decodeTracklets){
		return paths;
	}
	else{
		std::vector<std::vector<size_t>> origPaths(paths.size());
		for (int i = 0; i < paths.size(); ++i) {
			for (int j = 0; j < paths[i].size(); ++j) {
				size_t v=paths[i][j];
				for (int k = 0; k < smallPaths[v].size(); ++k) {
					size_t vertex=smallPaths[v][k];
					origPaths[i].push_back(vertex);
				}
			}

		}


		return origPaths;
	}
}


template<class T>
inline std::pair<double,double> Data<T>::evaluate(std::vector<std::vector<size_t>>& paths){
	//std::string infoFileName=parameters.getOutputFileName()+"-info.txt";
	//std::ofstream infoFile;
	//infoFile.open(infoFileName);

	readCompleteGraph();

	std::vector<size_t> labels(pGraphComplete->numberOfVertices(),0);
	for (int i = 0; i < paths.size(); ++i) {
		for (int j = 0; j < paths[i].size(); ++j) {
			size_t v=paths[i][j];
			labels[v]=i+1;
		}

	}

	double solverObjective=0;
	for (int i = 0; i < paths.size(); ++i) {
		for (int j = 0; j < paths[i].size()-1; ++j) {
			size_t v=paths[i][j];
			size_t w=paths[i][j+1];
			solverObjective+=(*pCompleteScore)[pGraphComplete->findEdge(v,w).second];
		}
	}



	double objValue=0;

	//std::map<size_t,std::map<size_t,double>> originalEdges;

	std::cout<<"Getting base edges from complete graph. "<<std::endl;
	parameters.infoFile()<<"Getting base edges from complete graph. "<<std::endl;
	parameters.infoFile().flush();
	//parameters.paramsFile<<"Getting base edges from complete graph. "<<std::endl;
	size_t maxLabel=0;
	for (int i = 0; i < pGraphComplete->numberOfEdges(); ++i) {

		size_t v=pGraphComplete->vertexOfEdge(i,0);
		size_t w=pGraphComplete->vertexOfEdge(i,1);
		//std::cout<<"edge "<<v<<" "<<w<<std::endl;
		size_t l0=pVertexGroups->getGroupIndex(v);
		size_t l1=pVertexGroups->getGroupIndex(w);
		maxLabel=std::max(maxLabel,labels[v]);
		maxLabel=std::max(maxLabel,labels[w]);

		if(labels[v]!=0&&labels[v]==labels[w]){
			double score=(*pCompleteScore)[i];
			if(score<0||l1-l0<=parameters.getRepulsiveTimeGap()){
				objValue+=score;
				if(l1-l0<=parameters.getTrTimeGapLifted()){
					solverObjective+=score;
				}
			}
		}

	}



	std::cout<<"clustering_evaluation "<<objValue<<std::endl;
	parameters.infoFile()<<"clustering_evaluation "<<objValue<<std::endl;

	//solverObjective+=objValue;
	solverObjective+=(parameters.getInputCost()+parameters.getOutputCost())*paths.size();
	std::cout<<"solver objective "<<solverObjective<<std::endl;
	parameters.infoFile()<<"solver objective "<<solverObjective<<std::endl;
	parameters.infoFile().flush();

	std::pair<double,double> toReturn(objValue,solverObjective);

	return toReturn;

}







template<class T>
inline std::unordered_map<size_t,std::vector<size_t>> Data<T>::findTrajectoryBreaks(std::vector<std::vector<size_t>>& paths){


	evaluate(paths);


	double inOutCost=parameters.getInputCost()+parameters.getOutputCost();
	std::unordered_map<size_t,std::vector<size_t>> outputBrakes;

		std::cout<<"finding track break points, number of tracks "<<paths.size()<<std::endl;
		for (int i = 0; i < paths.size(); ++i) {

			std::map<size_t,double> origBrakes;
			for (int j = 0; j < paths[i].size()-1; ++j) {
				size_t v=paths[i][j];
				size_t w=paths[i][j+1];
				auto fe=pGraphComplete->findEdge(v,w);
				origBrakes[j]=(*pCompleteScore)[fe.second];

			}

			for (int j = 0; j < paths[i].size()-1; ++j) {
				size_t v=paths[i][j];
				for (int k = j+1; k < paths[i].size(); ++k) {
					size_t w=paths[i][k];

					auto fe=pGraphComplete->findEdge(v,w);
					if(fe.first){
						double costToAdd=(*pCompleteScore)[fe.second];
						size_t timeW=pVertexGroups->getGroupIndex(w);
						size_t timeV=pVertexGroups->getGroupIndex(v);
						if(costToAdd<0||timeW-timeV<=parameters.getRepulsiveTimeGap()){
							if(timeW-timeV<=parameters.getTrTimeGapLifted()){
								for (int l = j; l < k; ++l) {
									origBrakes[l]+=costToAdd;
								}
							}
						}

					}
				}
			}

			for(auto tB:origBrakes){
				if(tB.second>inOutCost){
					outputBrakes[i].push_back(tB.first);
					std::cout<<"Not good track "<<i<<" better cut in vertex "<<paths[i][tB.first]<<", cost"<<tB.second<<std::endl;
					parameters.infoFile()<<"Not good track "<<i<<" better cut in vertex "<<paths[i][tB.first]<<", cost"<<tB.second<<std::endl;
				}
			}
			//std::cout<<std::endl;

		}


	std::cout<<"trajectory break points checked "<<std::endl;
	parameters.infoFile()<<"trajectory break points checked "<<std::endl;
	parameters.infoFile().flush();
	//parameters.paramsFile<<"time break points checked "<<std::endl;
	return outputBrakes;


}



template<class T>
inline std::unordered_map<size_t,std::set<size_t>> Data<T>::findTimeBreaks(std::vector<std::vector<size_t>>& paths,bool origOnly){


	evaluate(paths);

	double inOutCost=parameters.getInputCost()+parameters.getOutputCost();
	std::unordered_map<size_t,std::set<size_t>> outputBrakes;


	origOnly=true;

	if(!origOnly){//TODO make this general useful for checking if there is not an error in the solver, instead of timebreaks, return vertex breaks

		for (int i = 0; i < paths.size(); ++i) {
			std::map<size_t,double> timeBrakes;
			std::map<size_t,double> origTimeBrakes;
			for (int j = 0; j < paths[i].size()-1; ++j) {
				size_t v=paths[i][j];
				size_t w=paths[i][j+1];
				size_t time=pVertexGroups->getGroupIndex(v);
				auto findEdge=pGraph->findEdge(v,w);
				//			if(!findEdge.first){
				//				throw std::runtime_error(std::string("wrong connection"));
				//			}
				//			else{
				timeBrakes[time]=costs[getEdgeVarIndex(findEdge.second)];
				//origTimeBrakes[time]=costs[getEdgeVarIndex(findEdge.second)];
				origTimeBrakes[time]=0;
				//			}
			}

			for (int j = 0; j < paths[i].size()-1; ++j) {
				size_t v=paths[i][j];
				size_t timeV=pVertexGroups->getGroupIndex(v);
				for (int k = j+1; k < paths[i].size(); ++k) {
					size_t w=paths[i][k];
					auto findLiftedEdge=pGraphLifted->findEdge(v,w);
					size_t timeW=pVertexGroups->getGroupIndex(w);
					if(findLiftedEdge.first){
						size_t edge=findLiftedEdge.second;
						double edgeCost=costs[getLiftedEdgeVarIndex(edge)];
							for (int t = timeV; t < timeW; ++t) {
								if(timeBrakes.count(t)>0){
									timeBrakes[t]+=edgeCost;

								}
							}


					}
					auto fe=pGraphComplete->findEdge(v,w);
					if(fe.first){
						for (int t = timeV; t < timeW; ++t) {
							if(origTimeBrakes.count(t)>0){
								origTimeBrakes[t]+=(*pCompleteScore)[fe.second];

							}
						}

					}
				}
			}

			//std::cout<<"check track "<<i<<std::endl;
			for(auto tB:timeBrakes){
				if(tB.second>inOutCost){
					std::cout<<"Wrong track "<<i<<" better cut in time "<<tB.first<<"cost"<<tB.second<<std::endl;
					parameters.infoFile()<<"Wrong track "<<i<<" better cut in time "<<tB.first<<" cost "<<tB.second<<std::endl;
				}
			}
			for(auto tB:origTimeBrakes){
				if(tB.second>inOutCost){
					outputBrakes[i].insert(tB.first);
					std::cout<<"Not good track "<<i<<" better cut in time "<<tB.first<<" cost "<<tB.second<<std::endl;
					parameters.infoFile()<<"Not good track "<<i<<" better cut in time "<<tB.first<<" cost "<<tB.second<<std::endl;
				}
			}
			//std::cout<<std::endl;

		}

	}
	else{
		//TODO another value where direct neighbors are count two times (due to base edge),
		//the new value should be used to judge if a new computation makes sense
		//TODO maybe store just index in path instead of time
		std::cout<<"finding track break points, number of tracks "<<paths.size()<<std::endl;
		for (int i = 0; i < paths.size(); ++i) {

			std::map<size_t,double> origTimeBrakes;
			for (int j = 0; j < paths[i].size()-1; ++j) {
				size_t v=paths[i][j];
				size_t w=paths[i][j+1];
				size_t time=pVertexGroups->getGroupIndex(v);
				auto fe=pGraphComplete->findEdge(v,w);
				origTimeBrakes[time]=(*pCompleteScore)[fe.second];
//				origTimeBrakes[time]=0;

			}
			//if(paths[i].size()>1)	std::cout<<"path "<<(i+1)<<" has size "<<paths[i].size()<<std::endl;

			for (int j = 0; j < paths[i].size()-1; ++j) {
				size_t v=paths[i][j];
				size_t timeV=pVertexGroups->getGroupIndex(v);
				for (int k = j+1; k < paths[i].size(); ++k) {
					size_t w=paths[i][k];
					size_t timeW=pVertexGroups->getGroupIndex(w);

					auto fe=pGraphComplete->findEdge(v,w);
					if(fe.first){
						//double edgeCost=(*pCompleteScore)[fe.second];
						//if(edgeCost<=parameters.getNegativeThresholdLifted()||edgeCost>=parameters.getPositiveThresholdLifted()){
							for (int t = timeV; t < timeW; ++t) {
								if(origTimeBrakes.count(t)>0){
									double costToAdd=(*pCompleteScore)[fe.second];
									if(costToAdd<0||timeW-timeV<=parameters.getRepulsiveTimeGap()){
										if(timeW-timeV<=parameters.getTrTimeGapLifted()){
											origTimeBrakes[t]+=costToAdd;
										}
									}

								}
							}
						//}

					}
				}
			}

			for(auto tB:origTimeBrakes){
				if(tB.second>inOutCost){
					outputBrakes[i].insert(tB.first);
					std::cout<<"Not good track "<<i<<" better cut in time "<<tB.first<<"cost"<<tB.second<<std::endl;
					parameters.infoFile()<<"Not good track "<<i<<" better cut in time "<<tB.first<<"cost"<<tB.second<<std::endl;
				}
			}
			//std::cout<<std::endl;

		}

	}
	std::cout<<"time break points checked "<<std::endl;
	parameters.infoFile()<<"time break points checked "<<std::endl;
	parameters.infoFile().flush();
	//parameters.paramsFile<<"time break points checked "<<std::endl;
	return outputBrakes;


}




template<class T>
inline void Data<T>::outputSolution(std::vector<std::vector<size_t>>& paths,bool isIntervals){


	std::ofstream file;
	if(isIntervals){
		file.open(parameters.getOutputFileName() + "-all-paths-INTERVALS.txt");
	}
	else{
		file.open(parameters.getOutputFileName() + "-all-paths-FINAL.txt");
	}


	for (int i = 0; i < paths.size(); ++i) {
		//std::cout<<"output path "<<i<<std::endl;
		for (int j = 0; j < paths[i].size(); ++j) {
			size_t v=paths[i][j];
			file<<v<<" ";
		//	labels[v]=i+1;
		}
		file<<std::endl;
	}


	file.close();
	std::cout<<"file closed "<<std::endl;
	parameters.infoFile()<<"file closed "<<std::endl;
	parameters.infoFile().flush();




}




}

#endif /* INCLUDE_DISJOINT_PATHS_DISJOINT_PATHS_INIT_HXX_ */
