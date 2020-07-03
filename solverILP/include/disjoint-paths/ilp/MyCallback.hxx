/*
 * solver-disjoint-gurobi.hxx
 *
 *  Created on: Oct 16, 2019
 *      Author: fuksova
 */

#ifndef INCLUDE_DISJOINT_PATHS_ILP_MYCALLBACK_HXX_
#define INCLUDE_DISJOINT_PATHS_ILP_MYCALLBACK_HXX_


#include <cmath>
#include <vector>
#include <stack>
#include <sstream>
#include <fstream>
#include <iostream>
#include <map>
#include <list>
#include "disjoint-paths/disjointPathsData.hxx"
#include "disjoint-paths/disjointStructure.hxx"
#include "disjoint-paths/ilp/gurobi-callback-disjoint.hxx"





namespace disjointPaths {
namespace ilp {



class MyCallback: public disjointPaths::ilp::Gurobi::Callback
{
public:
	MyCallback(disjointPaths::ilp::Gurobi& solver, Data<>& data,bool outputs);
//	MyCallback(disjointPaths::ilp::Gurobi& solver, Data<>& data,bool experimental) :
//		disjointPaths::ilp::Gurobi::Callback(solver),
//		data_(data),
//		coefficients_(data.getCosts().size()),
//		variables_(data.getCosts().size()),
//		useExperimental(experimental)
//	{
//
//	}
	//~MyCallback() override;
	//~MyCallback() override{}


	void separateAndAddLazyConstraints() override;
	void simplePathConstraint(size_t w,size_t& counter,std::list<size_t> & path, std::vector<bool> &isOnPath,andres::graph::Digraph<> const & graph,andres::graph::Digraph<> const & liftedGraph);
	void simpleReversedPathConstraint(size_t w,size_t& counter,std::list<size_t> & path, std::vector<bool> &isOnPath,andres::graph::Digraph<> const & graph,andres::graph::Digraph<> const & liftedGraph);
	void simplePathToCutConstraint(size_t w,size_t& counter,std::list<size_t> & path, std::vector<bool> &isOnPath,andres::graph::Digraph<> const & graph,andres::graph::Digraph<> const & liftedGraph);
	void simpleCutConstraint(size_t w,size_t& counter,std::list<size_t> & path, std::vector<bool> &isOnPath,andres::graph::Digraph<> const & graph,andres::graph::Digraph<> const & liftedGraph);
	size_t preparePathConstraint(size_t vCurrent,size_t vertex,size_t& counter,std::unordered_map<size_t,size_t>& pred,std::unordered_map<size_t,bool>& predIsLifted, andres::graph::Digraph<> const & graph,andres::graph::Digraph<> const & liftedGraph,size_t w);
	void preparePathConstrReversed(size_t vCurrent,size_t vertex,size_t& counter,std::unordered_map<size_t,size_t>& pred,std::unordered_map<size_t,bool>& predIsLifted, andres::graph::Digraph<> const & graph,andres::graph::Digraph<> const & liftedGraph);
	void prepareCutReversed(size_t vCurrent,size_t vertex,size_t w,size_t& counter,std::unordered_map<size_t,size_t>& pred,std::unordered_map<size_t,bool>& predIsLifted, andres::graph::Digraph<> const & graph,andres::graph::Digraph<> const & liftedGraph);
	size_t addSimpleConstraints(bool useCutVersion);
	size_t addMixedCutConstraints();
	size_t addMixedConstraints();

protected:


	std::vector<double> coefficients_;
	Data<>& data_;
	size_t numberOfFeasibleSolutions_ { 0 };
	size_t numberOfSeparationCalls_ { 0 };
	std::vector<size_t> variables_;
	size_t allConstraints_{0};
	bool useOutputs;

};


//
//
//MyCallback::MyCallback(disjointPaths::ilp::Gurobi& solver, Data<>& data,bool experimental) :
//		disjointPaths::ilp::Gurobi::Callback(solver),
//		data_(data),
//		coefficients_(data.getCosts().size()),
//		variables_(data.getCosts().size()),
//		useExperimental(experimental)
//	{
//
//	}
//
//MyCallback::~MyCallback(){}
//
//void MyCallback::separateAndAddLazyConstraints(){
//
//	std::stringstream stream;
//
//	std::cout<<"separate and add lazy"<<std::endl;
//	data_.timer.stop();
//	auto time = data_.timer.get_elapsed_seconds();
//	data_.timer.start();
//
//	auto gap = (this->objectiveBest_ - this->objectiveBound_) / (1.0 + std::fabs(this->objectiveBest_));
//	std::cout << time
//			<< ' ' << this->objectiveBound_
//			<< ' ' << this->objectiveBest_
//			<< ' ' << gap
//			<< std::flush;
//	stream << time
//			<< ' ' << this->objectiveBound_
//			<< ' ' << this->objectiveBest_
//			<< ' ' << gap;
//
//	levinkov::Timer t_separation;
//	t_separation.start();
//
//
//
//	size_t numberOfAddedConstraints=0;
//
//	numberOfAddedConstraints=addMixedConstraints();
//	//numberOfAddedConstraints=addMixedCutConstraints();
//	//numberOfAddedConstraints=addSimpleConstraints(true);
//
//
//
//	std::cout << "added constraints " << numberOfAddedConstraints << std::flush;
//	stream << "added constraints " << numberOfAddedConstraints<<", ";
//	allConstraints_+=numberOfAddedConstraints;
//	stream << "all constraints " << allConstraints_<<", ";
//
//	std::cout <<std::endl;
//	stream <<std::endl;
//
//
//	t_separation.stop();
//	data_.timer.stop(); // not keeping time for writing log
//
//	auto objValue = .0;
//	for (size_t i = 0; i < data_.getCosts().size(); ++i){
//		objValue += data_.getCosts()[i]*this->label(i);
//	}
//
//
//	std::cout<<std::endl;
//
//	std::cout << "objective value " << objValue;
//	stream << ", " << objValue;
//
//	stream << ", " << t_separation.get_elapsed_seconds() << std::endl;
//	std::cout << "time " << t_separation.get_elapsed_seconds() << std::endl;
//
//	{
//		std::ofstream file(data_.parameters.getOutputFileName() + "-optimization-log.txt", std::ofstream::out | std::ofstream::app);
//		file << stream.str();
//		file.close();
//	}
//
//
//
//
//	if (numberOfAddedConstraints == 0)
//	{
//
//
//
//		std::ofstream file;
//
//		file.open(data_.parameters.getOutputFileName() + "-all-paths-FEASIBLE-" + std::to_string(numberOfFeasibleSolutions_) + ".txt");
//		for (int i = 0; i < data_.getGraph().numberOfEdgesFromVertex(data_.getSourceNode()); ++i) {
//			size_t edgeIndex=data_.getGraph().edgeFromVertex(data_.getSourceNode(),i);
//			if(this->label(data_.getEdgeVarIndex(edgeIndex))>0.5){
//				size_t v=data_.getGraph().vertexFromVertex(data_.getSourceNode(),i);
//				while(v!=data_.getTerminalNode()){
//					file<<v<<" ";
//					for (int j = 0; j < data_.getGraph().numberOfEdgesFromVertex(v); ++j) {
//						size_t edge=data_.getGraph().edgeFromVertex(v,j);
//						if(this->label(data_.getEdgeVarIndex(edge))>0.5){
//							v=data_.getGraph().vertexFromVertex(v,j);
//							break;
//						}
//					}
//				}
//				file<<std::endl;
//			}
//		}
//		file.close();
//
//
//
//		++numberOfFeasibleSolutions_;
//	}
//
//
//	std::cout<<"All lazy constraints: "<<allConstraints_<<std::endl;
//	++numberOfSeparationCalls_;
//
//	data_.timer.start(); // resume keeping time
//}
//
//
//
//void MyCallback::simplePathConstraint(size_t w,size_t& counter,std::list<size_t> & path, std::vector<bool> &isOnPath,andres::graph::Digraph<> const & graph,andres::graph::Digraph<> const & liftedGraph){
//		auto it=path.begin();
//		size_t v=*it;
//
//		for (int i = 0; i < graph.numberOfEdgesFromVertex(v); ++i) {
//			size_t u=graph.vertexFromVertex(v,i);
//			if(isOnPath[u]&&(data_.isReachable(u,w)||u==w)){
//				size_t e=graph.edgeFromVertex(v,i);
//				variables_[counter]=data_.getEdgeVarIndex(e);
//				coefficients_[counter]=1;
//				counter++;
//			}
//		}
//		it++;
//
//		while(*it!=w){
//			size_t vertexVarIndex=data_.getVertexVarIndex(*it);
//			variables_[counter]=data_.getVertexVarIndex(*it);
//			coefficients_[counter]=-1;
//			counter++;
//			for (int i = 0; i < graph.numberOfEdgesFromVertex(*it); ++i) {
//				size_t u=graph.vertexFromVertex(*it,i);
//				if(isOnPath[u]&&(data_.isReachable(u,w)||u==w)){
//					size_t e=graph.edgeFromVertex(*it,i);
//					variables_[counter]=data_.getEdgeVarIndex(e);
//					coefficients_[counter]=1;
//					counter++;
//				}
//			}
//			it++;
//		}
//	}
//
//
//
//	void MyCallback::simpleReversedPathConstraint(size_t w,size_t& counter,std::list<size_t> & path, std::vector<bool> &isOnPath,andres::graph::Digraph<> const & graph,andres::graph::Digraph<> const & liftedGraph){
//			auto it=path.begin();
//			size_t v=*it;
//
//			for (int i = 0; i < graph.numberOfEdgesFromVertex(v); ++i) {
//				size_t u=graph.vertexFromVertex(v,i);
//				if(isOnPath[u]&&(data_.isReachable(u,w)||u==w)){
//					size_t e=graph.edgeFromVertex(v,i);
//					variables_[counter]=data_.getEdgeVarIndex(e);
//					coefficients_[counter]=1;
//					counter++;
//				}
//			}
//			it++;
//
//			while(*it!=w){
//				for (int i = 0; i < graph.numberOfEdgesFromVertex(*it); ++i) {
//					size_t u=graph.vertexFromVertex(*it,i);
//					if(!isOnPath[u]||(!data_.isReachable(u,w)&&u!=w)){
//						size_t e=graph.edgeFromVertex(*it,i);
//						variables_[counter]=data_.getEdgeVarIndex(e);
//						coefficients_[counter]=-1;
//						counter++;
//					}
//				}
//				it++;
//			}
//		}
//
//
//
//	void MyCallback::simplePathToCutConstraint(size_t w,size_t& counter,std::list<size_t> & path, std::vector<bool> &isOnPath,andres::graph::Digraph<> const & graph,andres::graph::Digraph<> const & liftedGraph){
//		auto it=path.begin();
//		while(data_.isReachable(*it,w)){
//			size_t vertexVarIndex=data_.getVertexVarIndex(*it);
//			variables_[counter]=vertexVarIndex;
//			coefficients_[counter]=1;
//			counter++;
//			for (int i = 0; i < graph.numberOfEdgesFromVertex(*it); ++i) {
//				size_t u=graph.vertexFromVertex(*it,i);
//				if((isOnPath[u]&&data_.isReachable(u,w))||(!data_.isReachable(u,w)&&u!=w)){
//					size_t e=graph.edgeFromVertex(*it,i);
//					variables_[counter]=data_.getEdgeVarIndex(e);
//					coefficients_[counter]=-1;
//					counter++;
//				}
//			}
//			it++;
//		}
//	}
//
//
//
//
//	void MyCallback::simpleCutConstraint(size_t w,size_t& counter,std::list<size_t> & path, std::vector<bool> &isOnPath,andres::graph::Digraph<> const & graph,andres::graph::Digraph<> const & liftedGraph){
//		auto it=path.begin();
//		while(data_.isReachable(*it,w)){
//			for (int i = 0; i < graph.numberOfEdgesFromVertex(*it); ++i) {
//				size_t u=graph.vertexFromVertex(*it,i);
//				if(!isOnPath[u]&&(data_.isReachable(u,w)||u==w)){
//					size_t e=graph.edgeFromVertex(*it,i);
//					variables_[counter]=data_.getEdgeVarIndex(e);
//					coefficients_[counter]=1;
//					counter++;
//				}
//			}
//			it++;
//		}
//	}
//
//
//
//
//
//	size_t MyCallback::preparePathConstraint(size_t vCurrent,size_t vertex,size_t& counter,std::map<size_t,size_t>& pred,std::map<size_t,bool>& predIsLifted, andres::graph::Digraph<> const & graph,andres::graph::Digraph<> const & liftedGraph,size_t w){
//		std::vector<size_t> usedVertices;
//
//		size_t endVertex=vCurrent;
//		if(w!=vCurrent){
//			while(vCurrent!=vertex){
//				//std::cout<<vCurrent<<" ";
//
//				size_t newVertex=pred[vCurrent];
//				for (int i = 0; i < graph.numberOfEdgesFromVertex(newVertex); ++i) {
//
//					size_t v2=graph.vertexFromVertex(newVertex,i);
//					if(!data_.isReachable(v2,w)&&v2!=w){
//						size_t edgeIndex=graph.edgeFromVertex(newVertex,i);
//						size_t edgeVarIndex=data_.getEdgeVarIndex(edgeIndex);
//						variables_[counter]=edgeVarIndex;
//						coefficients_[counter]=-1;
//						counter++;
//					}
//				}
//
//
//				vCurrent=newVertex;
//			}
//
//		}
//		vCurrent=endVertex;
//		while(vCurrent!=vertex){
//			//std::cout<<vCurrent<<" ";
//			usedVertices.push_back(vCurrent);
//			size_t newVertex=pred[vCurrent];
//			size_t edgeIndex=0;
//			size_t edgeVarIndex=0;
//			if(predIsLifted[vCurrent]){
//
//				edgeIndex=liftedGraph.findEdge(newVertex,vCurrent).second;
//				edgeVarIndex=data_.getLiftedEdgeVarIndex(edgeIndex);
//
//			}
//			else{
//				edgeIndex=graph.findEdge(newVertex,vCurrent).second;
//				edgeVarIndex=data_.getEdgeVarIndex(edgeIndex);
//			}
//			variables_[counter]=edgeVarIndex;
//			coefficients_[counter]=-1;
//			counter++;
//			variables_[counter]=data_.getVertexVarIndex(newVertex);
//			coefficients_[counter]=1;
//			counter++;
//
//
//			vCurrent=newVertex;
//		}
//		size_t v=counter-1;
//		bool useCutLike=true;
//		if(useCutLike){
//			usedVertices.push_back(vertex);
//			for (int i = usedVertices.size()-1; i > 1; --i) {
//				//if(!predIsLifted[usedVertices[i-1]]){
//					for (int j = i-2; j >= 0; --j) {
//						auto findEdge=graph.findEdge(usedVertices[i],usedVertices[j]);
//						if(findEdge.first){
//							variables_[counter]=data_.getEdgeVarIndex(findEdge.second);
//							coefficients_[counter]=-1;
//							counter++;
//						}
//					}
//
//				//}
//			}
//		}
//
//
//		return v;
//	}
//
//	//TODO implement cut version of inequalities using P'
//	void MyCallback::preparePathConstrReversed(size_t vCurrent,size_t vertex,size_t& counter,std::map<size_t,size_t>& pred,std::map<size_t,bool>& predIsLifted, andres::graph::Digraph<> const & graph,andres::graph::Digraph<> const & liftedGraph){
//			//std::set<size_t> usedVertices;
//			std::vector<bool> used(graph.numberOfVertices(),0);
//			size_t edgeIndex=0;
//			size_t edgeVarIndex=0;
//
//			while(vCurrent!=vertex){
//				//std::cout<<vCurrent<<" ";
//				//usedVertices.insert(vCurrent);
//				used[vCurrent]=1;
//				size_t newVertex=pred[vCurrent];
//				if(newVertex==vertex) break;
//
//				if(predIsLifted[vCurrent]){
//					edgeIndex=liftedGraph.findEdge(newVertex,vCurrent).second;
//					edgeVarIndex=data_.getLiftedEdgeVarIndex(edgeIndex);
//					variables_[counter]=edgeVarIndex;
//					coefficients_[counter]=1;
//					counter++;
//					auto findEdge=graph.findEdge(newVertex,vCurrent);
//					if(findEdge.first){
//						edgeIndex=findEdge.second;
//						edgeVarIndex=data_.getEdgeVarIndex(edgeIndex);
//						variables_[counter]=edgeVarIndex;
//						coefficients_[counter]=-1;
//						counter++;
//					}
//				}
//				for (int i = 0; i < graph.numberOfEdgesFromVertex(newVertex); ++i) {
//					size_t u=graph.vertexFromVertex(newVertex,i);
//					//if(usedVertices.count(u)==0){
//					if(!used[u]){
//						edgeIndex=graph.edgeFromVertex(newVertex,i);
//						edgeVarIndex=data_.getEdgeVarIndex(edgeIndex);
//						variables_[counter]=edgeVarIndex;
//						coefficients_[counter]=-1;
//						counter++;
//
//					}
//
//				}
//				vCurrent=newVertex;
//			}
//			bool firstIsLifted=predIsLifted[vCurrent];
//			for (int i = 0; i < graph.numberOfEdgesFromVertex(vertex); ++i) {
//				size_t u=graph.vertexFromVertex(vertex,i);
//				//if(usedVertices.count(u)>0&&(u!=vCurrent||!firstIsLifted)){
//				if(used[u]&&(u!=vCurrent||!firstIsLifted)){
//					edgeIndex=graph.edgeFromVertex(vertex,i);
//					edgeVarIndex=data_.getEdgeVarIndex(edgeIndex);
//					variables_[counter]=edgeVarIndex;
//					coefficients_[counter]=1;
//					counter++;
//
//				}
//
//			}
//			if(firstIsLifted){
//				edgeIndex=liftedGraph.findEdge(vertex,vCurrent).second;
//				edgeVarIndex=data_.getLiftedEdgeVarIndex(edgeIndex);
//				variables_[counter]=edgeVarIndex;
//				coefficients_[counter]=1;
//				counter++;
//			}
//
//		}
//
//
//
//	void MyCallback::prepareCutReversed(size_t vCurrent,size_t vertex,size_t w,size_t& counter,std::map<size_t,size_t>& pred,std::map<size_t,bool>& predIsLifted, andres::graph::Digraph<> const & graph,andres::graph::Digraph<> const & liftedGraph){
//		//std::set<size_t> usedVertices;
//		std::vector<bool> used(graph.numberOfVertices(),0);
//
//		while(vCurrent!=vertex){
//			//std::cout<<vCurrent<<" ";
//			//usedVertices.insert(vCurrent);
//			used[vCurrent]=1;
//
//			size_t newVertex=pred[vCurrent];
//			size_t edgeIndex=0;
//			size_t edgeVarIndex=0;
//			if(predIsLifted[vCurrent]){
//				edgeIndex=liftedGraph.findEdge(newVertex,vCurrent).second;
//				edgeVarIndex=data_.getLiftedEdgeVarIndex(edgeIndex);
//				variables_[counter]=edgeVarIndex;
//				coefficients_[counter]=-1;
//				counter++;
//				auto findEdge=graph.findEdge(newVertex,vCurrent);
//				if(findEdge.first){
//					edgeIndex=findEdge.second;
//					edgeVarIndex=data_.getEdgeVarIndex(edgeIndex);
//					variables_[counter]=edgeVarIndex;
//					coefficients_[counter]=1;
//					counter++;
//				}
//			}
//			for (int i = 0; i < graph.numberOfEdgesFromVertex(newVertex); ++i) {
//				size_t u=graph.vertexFromVertex(newVertex,i);
//				//if(usedVertices.count(u)==0&&(data_.isReachable(u,w)||u==w)){
//				if(!used[u]&&(data_.isReachable(u,w)||u==w)){
//					edgeIndex=graph.edgeFromVertex(newVertex,i);
//					edgeVarIndex=data_.getEdgeVarIndex(edgeIndex);
//					variables_[counter]=edgeVarIndex;
//					coefficients_[counter]=1;
//					counter++;
//
//				}
//
//			}
//			vCurrent=newVertex;
//		}
//	}
//
//
//
//	size_t MyCallback::addSimpleConstraints(bool useCutVersion){
//		size_t constraintsCounter=0;
//		size_t t=data_.getTerminalNode();
//		size_t s=data_.getSourceNode();
//		andres::graph::Digraph<> const& graph=data_.getGraph();  //TODO check if the const & is correct
//		andres::graph::Digraph<> const& liftedGraph=data_.getGraphLifted();
//
//		//bool useCutVersion=false;
//
//
//		for (int i = 0; i < graph.numberOfEdgesToVertex(t); ++i) {
//
//
//			size_t edge=graph.edgeToVertex(t,i);
//			size_t eVarIndex=data_.getEdgeVarIndex(edge);
//			if(this->label(eVarIndex) > 0.5){
//				std::vector<bool> isOnPath(graph.numberOfVertices(),0);
//				size_t vertex=graph.vertexToVertex(t,i);
//				std::list<size_t> path;
//
//				path.push_front(t);
//				path.push_front(vertex);
//				isOnPath[vertex]=1;
//				isOnPath[t]=1;
//
//
//				size_t oldVertex=t;
//
//				while(vertex!=s){
//
//					for (int j = 0; j < liftedGraph.numberOfEdgesFromVertex(vertex); ++j) {
//						size_t liftedEdge=liftedGraph.edgeFromVertex(vertex,j);
//						size_t w=liftedGraph.vertexFromVertex(vertex,j);
//						size_t leVarIndex=data_.getLiftedEdgeVarIndex(liftedEdge);
//						bool isActive=this->label(leVarIndex)>0.5;
//						//std::cout<<" lifted "<<liftedEdge<<" active: "<<isActive<<std::endl;
//
//						if(w==oldVertex){
//
//							if(!isActive){
//								throw std::runtime_error(std::string("Missing lifted covering base constr. "));
//							}
//
//						}else if(!data_.isReachable(oldVertex,w)){
//
//							if(isActive){
//								throw std::runtime_error(std::string("Missing lifted unreachable base constr. "));
//							}
//
//						}else{//w is reachable from oldVertex
//
//
//							if(isOnPath[w]&&!isActive){
//
//								size_t counter=0;
//								size_t vCurrent=w;
//
//								if(useCutVersion){
//									simpleReversedPathConstraint(w,counter,path,isOnPath,graph,liftedGraph);
//								}
//								else{
//									simplePathConstraint(w,counter,path,isOnPath,graph,liftedGraph);
//								}
//
//								variables_[counter]=leVarIndex;
//								coefficients_[counter]=-1;
//
//								this->addLazyConstraint(variables_.begin(),	variables_.begin() + counter +1,
//										coefficients_.begin(),-std::numeric_limits<double>::infinity(),0);
//								constraintsCounter++;
//								//std::cout<<"positive lifted added "<<std::endl;
//
//							}
//
//							else if(!isOnPath[w]&&isActive){
//
//								size_t counter=0;
//								if(useCutVersion){
//									simpleCutConstraint(w,counter,path,isOnPath,graph,liftedGraph);
//								}
//								else{
//									simplePathToCutConstraint(w,counter,path,isOnPath,graph,liftedGraph);
//								}
//
//								variables_[counter]=leVarIndex;
//								coefficients_[counter]=-1;
//
//								this->addLazyConstraint(variables_.begin(),	variables_.begin() + (counter+1),
//										coefficients_.begin(),0,std::numeric_limits<double>::infinity());
//								constraintsCounter++;
//
//							}
//
//						}
//
//					}
//
//					size_t newVertex=0;
//					for (int j = 0; j < graph.numberOfEdgesToVertex(vertex);
//							++j) {
//						edge = graph.edgeToVertex(vertex, j);
//						eVarIndex = data_.getEdgeVarIndex(edge);
//
//						if (this->label(eVarIndex) > 0.5) {
//
//							newVertex = graph.vertexToVertex(vertex, j);
//							path.push_front(newVertex);
//							isOnPath[newVertex]=1;
//
//							break;
//						}
//					}
//					oldVertex=vertex;
//					vertex=newVertex;
//				}
//				while(!path.empty()){
//					isOnPath[path.front()]=0;
//					path.pop_front();
//				}
//
//			}
//		}
//
//
//		return constraintsCounter;
//	}
//
//
//
//	size_t MyCallback::addMixedCutConstraints(){
//		size_t constraintsCounter=0;
//		size_t t=data_.getTerminalNode();
//		size_t s=data_.getSourceNode();
//		andres::graph::Digraph<> const& graph=data_.getGraph();  //TODO check if the const & is correct
//		andres::graph::Digraph<> const& liftedGraph=data_.getGraphLifted();
//
//
//
//		for (int i = 0; i < graph.numberOfEdgesToVertex(t); ++i) {
//
//
//			size_t edge=graph.edgeToVertex(t,i);
//			size_t eVarIndex=data_.getEdgeVarIndex(edge);
//			if(this->label(eVarIndex) > 0.5){
//				std::vector<bool> isOnPath(graph.numberOfVertices(),0);
//				std::map<size_t,size_t> pred;
//				std::map<size_t,bool> predIsLifted;
//				std::map<size_t,size_t> predNeg;
//				size_t vertex=graph.vertexToVertex(t,i);
//				std::list<size_t> path;
//
//				path.push_front(t);
//				path.push_front(vertex);
//				isOnPath[vertex]=1;
//				isOnPath[t]=1;
//
//
//
//				size_t oldVertex=t;
//
//				while(vertex!=s){
//					pred[oldVertex]=vertex;
//					predIsLifted[oldVertex]=false;
//				//	std::cout<<"vertex "<<vertex<<std::endl;
//
//					for (int j = 0; j < liftedGraph.numberOfEdgesFromVertex(vertex); ++j) {
//						size_t liftedEdge=liftedGraph.edgeFromVertex(vertex,j);
//						size_t w=liftedGraph.vertexFromVertex(vertex,j);
//						size_t leVarIndex=data_.getLiftedEdgeVarIndex(liftedEdge);
//						bool isActive=this->label(leVarIndex)>0.5;
//						//std::cout<<" lifted "<<liftedEdge<<" active: "<<isActive<<std::endl;
//
//						if(w==oldVertex){
//
//							if(!isActive){
//								throw std::runtime_error(std::string("Missing lifted covering base constr. "));
//							}
//
//						}else if(!data_.isReachable(oldVertex,w)){
//
//							if(isActive){
//								throw std::runtime_error(std::string("Missing lifted unreachable base constr. "));
//							}
//
//						}else{//w is reachable from oldVertex
//
//
//								if(isOnPath[w]&&!isActive){
//
//									size_t counter=0;
//									size_t vCurrent=w;
//
//									preparePathConstrReversed(vCurrent,vertex,counter,pred,predIsLifted,graph,liftedGraph);
//
//									variables_[counter]=leVarIndex; //replace the last vertex variable with lifted edge variable
//									coefficients_[counter]=-1;
//
//									this->addLazyConstraint(variables_.begin(),	variables_.begin() + counter +1,
//											coefficients_.begin(),-std::numeric_limits<double>::infinity(),0);
//									constraintsCounter++;
//									//std::cout<<"positive lifted added "<<std::endl;
//
//								}
//
//							    else if(!isOnPath[w]&&isActive){
//									//
//									if(predNeg.count(w)>0){
//										size_t lastVertex=predNeg[w];
//										size_t counter=0;
//										size_t vCurrent=lastVertex;
//
//										prepareCutReversed(vCurrent,vertex,w,counter,pred,predIsLifted,graph,liftedGraph);
//										//std::cout<<std::endl;
//										size_t lastEdge=liftedGraph.findEdge(lastVertex,w).second;
//										size_t lastEdgeVarIndex=data_.getLiftedEdgeVarIndex(lastEdge);
//										variables_[counter]=lastEdgeVarIndex;
//										coefficients_[counter]=1;
//
//										counter++;
//
//										variables_[counter]=leVarIndex;
//										coefficients_[counter]=-1;
//										this->addLazyConstraint(variables_.begin(),	variables_.begin() + (counter+1),
//												coefficients_.begin(),0,std::numeric_limits<double>::infinity());
//										constraintsCounter++;
//
//
//									}
//									else{
//
//										auto itv=path.begin();
//										while(data_.isReachable(*itv,w)){
//											itv++;
//										}
//
//										size_t counter=0;
//										itv--;
//										size_t vCurrent=*itv;
//
//
//										prepareCutReversed(vCurrent,vertex,w,counter,pred,predIsLifted,graph,liftedGraph);
//
//										for (int i = 0; i < graph.numberOfEdgesFromVertex(vCurrent); ++i) {
//											size_t edge=graph.edgeFromVertex(vCurrent,i);
//											size_t v2=graph.vertexFromVertex(vCurrent,i);
//											if(data_.isReachable(v2,w)||v2==w){
//												size_t edgeVarIndex=data_.getEdgeVarIndex(edge);
//												variables_[counter]=edgeVarIndex;
//												coefficients_[counter]=1;
//												counter++;
//											}
//
//										}
//
//										variables_[counter]=leVarIndex;
//										coefficients_[counter]=-1;
//
//										this->addLazyConstraint(variables_.begin(),	variables_.begin() + (counter+1),
//												coefficients_.begin(),0,std::numeric_limits<double>::infinity());
//										constraintsCounter++;
//
//									}
//
//
//							    }
//
//
//						}
//						if(isOnPath[w]&&isActive){
//							pred[w]=vertex;
//							predIsLifted[w]=true;
//						}
//						if(!isOnPath[w]&&!isActive){
//							predNeg[w]=vertex;
//						}
//					}
//
//
//
//					size_t newVertex=0;
//					for (int j = 0; j < graph.numberOfEdgesToVertex(vertex);
//							++j) {
//						edge = graph.edgeToVertex(vertex, j);
//						eVarIndex = data_.getEdgeVarIndex(edge);
//
//						if (this->label(eVarIndex) > 0.5) {
//
//							newVertex = graph.vertexToVertex(vertex, j);
//							path.push_front(newVertex);
//							isOnPath[newVertex]=1;
//
//							break;
//						}
//					}
//					oldVertex=vertex;
//					vertex=newVertex;
//				}
//				while(!path.empty()){
//					isOnPath[path.front()]=0;
//					path.pop_front();
//				}
//				pred.clear();
//				predIsLifted.clear();
//				predNeg.clear();
//
//
//
//			}
//		}
//
//
//		return constraintsCounter;
//	}










}
}




#endif /* INCLUDE_DISJOINT_PATHS_ILP_MYCALLBACK_HXX_ */
