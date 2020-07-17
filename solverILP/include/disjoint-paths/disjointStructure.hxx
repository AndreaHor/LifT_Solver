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
#include <utility>

namespace disjointPaths {

template<class T = size_t>
struct DisjointStructure {
public:
	typedef T size_type;


	DisjointStructure(disjointPaths::DisjointParams<>& configParameters,char delim=',',CompleteStructure<>* cs=0,size_t minTime=0,size_t maxTime=0);
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


	disjointPaths::DisjointParams<>& parameters;

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
inline DisjointStructure<T>::DisjointStructure(disjointPaths::DisjointParams<>& configParameters,char delim,CompleteStructure<>* cs,size_t minTime,size_t maxTime):
parameters(configParameters)
{

	std::cout<<"interval "<<minTime<<","<<maxTime<<std::endl;
	parameters.infoFile()<<"interval "<<minTime<<","<<maxTime<<std::endl;
	parameters.infoFile().flush();

	useTimeFrames=parameters.isRestrictFrames()||parameters.isSparsify();
	size_t maxVertex;
	if(cs==0){
	//	if(useTimeFrames){
            vertexGroups=VertexGroups<size_t>(configParameters);
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


}





#endif /* INCLUDE_DISJOINT_PATHS_DISJOINTSTRUCTURE_HXX_ */
