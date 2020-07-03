/*
 * solver-disjoint-ilp.hxx
 *
 *  Created on: Sep 10, 2018
 *      Author: fuksova
 */

#ifndef INCLUDE_DISJOINT_PATHS_ILP_SOLVER_DISJOINT_ILP_HXX_
#define INCLUDE_DISJOINT_PATHS_ILP_SOLVER_DISJOINT_ILP_HXX_


#include <disjoint-paths/ilp/MyCallback.hxx>
#include <cmath>
#include <vector>
#include <stack>
#include <sstream>
#include <fstream>
#include <iostream>
#include <map>
#include <list>
#include "disjoint-paths/disjoint-paths-init.hxx"
#include "disjoint-paths/ilp/MyCallback.hxx"
#include "disjoint-paths/disjointStructure.hxx"
#include "disjoint-paths/disjointPathsMethods.hxx"





namespace disjointPaths {


template<class T=size_t>
void ilp_add_constraints(ilp::Gurobi& ilp,
		Data<>& data)
{


	std::vector<size_t> variables=std::vector<size_t>(data.numberOfEdges);
	std::vector<double> coeffs=std::vector<double>(data.numberOfEdges);
	andres::graph::Digraph<> const & graph=data.getGraph();
	andres::graph::Digraph<> const & graphLifted=data.getGraphLifted();

	for (int v = 0; v < data.numberOfVertices; ++v) {

		if(v!=data.getSourceNode()&&v!=data.getTerminalNode()){
			//all ineq. from the first paragraph except the path constraints
			size_t inflow = graph.numberOfEdgesToVertex(v);
			size_t outflow = graph.numberOfEdgesFromVertex(v);
			for (int j = 0; j < inflow; ++j) {

				size_t edge = graph.edgeToVertex(v,j);
				coeffs[j] = 1;
				variables[j] = data.getEdgeVarIndex(edge);

			}
			for (int j = 0; j < outflow; ++j) {
				size_t edge = graph.edgeFromVertex(v, j);
				coeffs[j+inflow] = -1;
				variables[j+inflow] = data.getEdgeVarIndex(edge);

			}

			//flow constraint for standard edges

			ilp.addConstraint(variables.begin(), variables.begin() + (inflow+outflow),
					coeffs.begin(), 0, 0);



			variables[inflow]=data.getVertexVarIndex(v);
			coeffs[inflow]=-1;


			//active vertex iff there is a flow
			ilp.addConstraint(variables.begin(), variables.begin() + (inflow+1),
					coeffs.begin(), 0, 0);


			//Constraints for lifted edges
			size_t numberOfLifted=graphLifted.numberOfEdgesFromVertex(v);


			//Without HOF, bad initial solution and slover computation
			bool useHigherOrderFlow=data.useTime();
			//bool useHigherOrderFlow=useHOF;

			std::map<size_t,std::vector<size_t>> groupToVertex;
			for (int i = 0; i < numberOfLifted; ++i) {
				size_t vertex2=graphLifted.vertexFromVertex(v,i);
				if(data.useTime()){
					size_t gi=data.getGroupIndex(vertex2);
					groupToVertex[gi].push_back(vertex2);
				}
				else if(data.useTimeForTracks()){
					const std::unordered_set<size_t>& allGroupsToV2=data.getTrackletsToGroups()[vertex2];
					for(size_t timeID:allGroupsToV2){
						groupToVertex[timeID].push_back(vertex2);
					}
				}
				//				size_t edge=data.getGraphLifted().edgeFromVertex(v,i);
				//				size_t leIndex=data.getLiftedEdgeVarIndex(edge);
				//				variables[0]=leIndex;
				//				coeffs[0]=-1;
				//				variables[1]=data.getVertexVarIndex(v);
				//				coeffs[1]=1;
				//
				//				ilp.addConstraint(variables.begin(), variables.begin() + 2,
				//						coeffs.begin(), 0, std::numeric_limits<double>::infinity());

			}

			if(data.useTime()||data.useTimeForTracks()){
				for(auto it=groupToVertex.begin();it!=groupToVertex.end();it++){
					std::vector<size_t>& groupVertices=it->second;
					if(groupVertices.size()>1){
						for (int i = 0; i < groupVertices.size(); ++i) {
							size_t le=graphLifted.findEdge(v,groupVertices[i]).second;
							variables[i]=data.getLiftedEdgeVarIndex(le);
							coeffs[i]=1;

						}
						variables[groupVertices.size()]=data.getVertexVarIndex(v);
						coeffs[groupVertices.size()]=-1;
						ilp.addConstraint(variables.begin(), variables.begin() + groupVertices.size()+1,
								coeffs.begin(), - std::numeric_limits<double>::infinity(),0);
					}

				}
				groupToVertex.clear();
			}
			for (int i = 0; i < graphLifted.numberOfEdgesToVertex(v); ++i) {
				size_t vertex2=graphLifted.vertexToVertex(v,i);
				if(data.useTime()){
					size_t gi=data.getGroupIndex(vertex2);
					groupToVertex[gi].push_back(vertex2);
				}
				else if(data.useTimeForTracks()){
					const std::unordered_set<size_t>& allGroupsToV2=data.getTrackletsToGroups()[vertex2];
					for(size_t timeID:allGroupsToV2){
						groupToVertex[timeID].push_back(vertex2);
					}
				}
				//				size_t edge=data.getGraphLifted().edgeToVertex(v,i);
				//				size_t leIndex=data.getLiftedEdgeVarIndex(edge);
				//				variables[0]=leIndex;
				//				coeffs[0]=-1;
				//				variables[1]=data.getVertexVarIndex(v);
				//				coeffs[1]=1;
				//
				//				ilp.addConstraint(variables.begin(), variables.begin() + 2,
				//						coeffs.begin(), 0, std::numeric_limits<double>::infinity());

			}
			if(data.useTime()||data.useTimeForTracks()){
				for(auto it=groupToVertex.begin();it!=groupToVertex.end();it++){
					std::vector<size_t>& groupVertices=it->second;
					if(groupVertices.size()>1){
						for (int i = 0; i < groupVertices.size(); ++i) {
							size_t le=graphLifted.findEdge(groupVertices[i],v).second;
							variables[i]=data.getLiftedEdgeVarIndex(le);
							coeffs[i]=1;

						}
						variables[groupVertices.size()]=data.getVertexVarIndex(v);
						coeffs[groupVertices.size()]=-1;
						ilp.addConstraint(variables.begin(), variables.begin() + groupVertices.size()+1,
								coeffs.begin(), - std::numeric_limits<double>::infinity(),0);
					}

				}
			}


			//}


		}
	}

	//Constraints for lifted covering base edge
	for (int e = 0; e < graph.numberOfEdges(); ++e) {
		size_t v0=graph.vertexOfEdge(e,0);
		size_t v1=graph.vertexOfEdge(e,1);

		auto findEdge=graphLifted.findEdge(v0,v1);
		if(findEdge.first){
			size_t leIndex=data.getLiftedEdgeVarIndex(findEdge.second);
			variables[0]=leIndex;
			coeffs[0]=1;
			variables[1]=data.getEdgeVarIndex(e);
			coeffs[1]=-1;
			ilp.addConstraint(variables.begin(), variables.begin() + 2,
					coeffs.begin(), 0, std::numeric_limits<double>::infinity());
		}

	}

	//Constraints for lifted unreachable from base
	for (int e = 0; e < graphLifted.numberOfEdges(); ++e) {
		size_t v0=graphLifted.vertexOfEdge(e,0);
		size_t v1=graphLifted.vertexOfEdge(e,1);
		size_t leIndex=data.getLiftedEdgeVarIndex(e);

		size_t counter=0;
		for (int i = 0; i < graph.numberOfEdgesFromVertex(v0); ++i) {
			size_t w=graph.vertexFromVertex(v0,i);
			if(w!=v1&&!data.isReachable(w,v1)){
				size_t edge=graph.edgeFromVertex(v0,i);
				size_t edgeIndex=data.getEdgeVarIndex(edge);
				variables[counter]=edgeIndex;
				coeffs[counter]=1;
				counter++;
			}

		}

		//if(counter>0){
		variables[counter]=leIndex;
		coeffs[counter]=1;
		counter++;
		variables[counter]=data.getVertexVarIndex(v0);
		coeffs[counter]=-1;
		ilp.addConstraint(variables.begin(), variables.begin() + counter+1,
				coeffs.begin(), -std::numeric_limits<double>::infinity(),0);
		//}
		counter=0;
		for (int i = 0; i < graph.numberOfEdgesToVertex(v1); ++i) {
			size_t w=graph.vertexToVertex(v1,i);
			if(w!=v0&&!data.isReachable(v0,w)){
				size_t edge=graph.edgeToVertex(v1,i);
				size_t edgeIndex=data.getEdgeVarIndex(edge);
				variables[counter]=edgeIndex;
				coeffs[counter]=1;
				counter++;
			}

		}
		//if(counter>0){
		variables[counter]=leIndex;
		coeffs[counter]=1;
		counter++;
		variables[counter]=data.getVertexVarIndex(v1);
		coeffs[counter]=-1;
		ilp.addConstraint(variables.begin(), variables.begin() + counter+1,
				coeffs.begin(), -std::numeric_limits<double>::infinity(),0);
		//}

	}





}


template<class T=size_t>
std::vector<double> ilp_init_solution(Data<>& data,bool outputs=true)
{
	disjointPaths::ilp::Gurobi ilp0;

	ilp0.setVerbosity(true);

	ilp0.setRelativeGap(0.0);
	ilp0.setAbsoluteGap(0.0);

	if(!outputs){
		ilp0.setOutputFlag(0);
		//ilp.setLogFile(data.parameters.getOutputFileName() +"-gurobi-temp.log");
	}


	ilp0.addVariables(data.getCosts().size(),data.getCosts().data());


	ilp_add_constraints<T>(ilp0,data);
	ilp0.optimize();

	if(outputs){
		std::cout<<"init solution"<<std::endl;
		data.parameters.infoFile()<<"init solution"<<std::endl;
		data.parameters.infoFile().flush();
	}
	std::vector<double> solution0(data.getCosts().size(),0);
	for (int i = 0; i < solution0.size(); ++i) {
		solution0[i]=ilp0.label(i);
		//if(ilp0.label(i)>0.5) std::cout<<ilp0.label(i)<<std::endl;
	}

	double costDifference=0;


	if(outputs){
		std::cout<<"fix init solution"<<std::endl;
		data.parameters.infoFile()<<"fix init solution"<<std::endl;
		data.parameters.infoFile().flush();
	}
	size_t t=data.getTerminalNode();
	size_t s=data.getSourceNode();
	andres::graph::Digraph<> const& graph=data.getGraph();  //TODO check if the const & is correct
	andres::graph::Digraph<> const& liftedGraph=data.getGraphLifted();
	std::vector<bool> isOnPath(graph.numberOfVertices(),0);

	//std::vector<size_t> pathToEnd(graph.numberOfVertices());
	size_t activeOutsidePath=0;
	size_t inactiveOnPath=0;

	for (int i = 0; i < graph.numberOfEdgesToVertex(t); ++i) {
		//std::cout<<i<<"-th edge to t"<<std::endl;

		size_t e=graph.edgeToVertex(t,i);
		size_t eVarIndex=data.getEdgeVarIndex(e);
		if(ilp0.label(eVarIndex) > 0.5){
			//std::map<size_t,size_t> vertexToDepth;
			size_t vertex=graph.vertexToVertex(t,i);
			isOnPath[vertex]=1;
			isOnPath[t]=1;  //Maybe not necessary

			std::list<size_t> path;
			path.push_front(t);
			path.push_front(vertex);

			while(vertex!=s){

				for (int j = 0; j < liftedGraph.numberOfEdgesFromVertex(vertex); ++j) {
					size_t le=liftedGraph.edgeFromVertex(vertex,j);
					size_t vertex2=liftedGraph.vertexFromVertex(vertex,j);
					size_t leVarIndex=data.getLiftedEdgeVarIndex(le);
					//size_t vertex2VarIndex=data_.getVertexVarIndex(vertex2);
					if(isOnPath[vertex2]) {
						solution0[leVarIndex]=1;
						if(ilp0.label(leVarIndex)==0){
							inactiveOnPath++;
							costDifference+=data.getCosts()[leVarIndex];
						}
					}
					else{
						solution0[leVarIndex]=0;
						if(ilp0.label(leVarIndex)==1){
							activeOutsidePath++;
							costDifference-=data.getCosts()[leVarIndex];
						}
					}


				}


				bool nextFound=false;
				size_t newVertex=0;

				for (int j = 0; j < graph.numberOfEdgesToVertex(vertex);
						++j) {
					e = graph.edgeToVertex(vertex, j);
					eVarIndex = data.getEdgeVarIndex(e);

					if (ilp0.label(eVarIndex) > 0.5) {
						if(nextFound){
							//							auto it=path.begin();
							//							size_t first=*it;
							//							size_t second=graph.vertexToVertex(vertex,j);
							//							it++;
							std::cout<<"Init: Two active incoming edges to vertex "<<std::endl;
							data.parameters.infoFile()<<"Init: Two active incoming edges to vertex "<<std::endl;
						}
						else{
							newVertex = graph.vertexToVertex(vertex, j);
							path.push_front(newVertex);
							isOnPath[newVertex]=1;

							nextFound=true;
						}

						//break;
					}
				}
				vertex=newVertex;



			}

			//std::cout<<std::endl;
			while(!path.empty()){
				isOnPath[path.front()]=0;
				path.pop_front();
			}
		}
	}
	if(outputs){
		std::cout<<"Init solution: inactive on path: "<<inactiveOnPath<<", active outside path: "<<activeOutsidePath<<std::endl;
		std::cout<<"cost difference "<<costDifference<<std::endl;
		data.parameters.infoFile()<<"Init solution: inactive on path: "<<inactiveOnPath<<", active outside path: "<<activeOutsidePath<<std::endl;
		data.parameters.infoFile()<<"cost difference "<<costDifference<<std::endl;
		data.parameters.infoFile().flush();
	}
	return solution0;



}



template<class T=size_t>
void ilp_check_solution(Data<>& data,std::vector<double>& solution0,bool outputs=true)
{

	if(outputs){
		std::cout<<"fix first solution"<<std::endl;
		data.parameters.infoFile()<<"fix first solution"<<std::endl;
		data.parameters.infoFile().flush();
	}
	size_t t=data.getTerminalNode();
	size_t s=data.getSourceNode();
	andres::graph::Digraph<> const& graph=data.getGraph();  //TODO check if the const & is correct
	andres::graph::Digraph<> const& liftedGraph=data.getGraphLifted();
	std::vector<bool> isOnPath(graph.numberOfVertices(),0);
	std::vector<bool> isOnAnyPath(graph.numberOfVertices(),0);

	//std::vector<size_t> pathToEnd(graph.numberOfVertices());

	size_t trackCounter=0;
	size_t activeEdgeCounter=0;
	for (int i = 0; i < graph.numberOfEdgesToVertex(t); ++i) {
		//std::cout<<i<<"-th edge to t"<<std::endl;

		size_t e=graph.edgeToVertex(t,i);
		size_t eVarIndex=data.getEdgeVarIndex(e);
		if(solution0[eVarIndex] > 0.5){
			//std::cout<<"check track "<<trackCounter<<std::endl;
			trackCounter++;
			activeEdgeCounter++;
			//std::map<size_t,size_t> vertexToDepth;
			size_t vertex=graph.vertexToVertex(t,i);
			isOnPath[vertex]=1;
			isOnPath[t]=1;  //Maybe not necessary
			isOnAnyPath[vertex]=1;
			isOnAnyPath[t]=1;

			std::list<size_t> path;
			path.push_front(t);
			path.push_front(vertex);

			while(vertex!=s){

				for (int j = 0; j < liftedGraph.numberOfEdgesFromVertex(vertex); ++j) {
					size_t le=liftedGraph.edgeFromVertex(vertex,j);
					size_t vertex2=liftedGraph.vertexFromVertex(vertex,j);
					size_t leVarIndex=data.getLiftedEdgeVarIndex(le);
					//size_t vertex2VarIndex=data_.getVertexVarIndex(vertex2);
					if(isOnPath[vertex2]&&solution0[leVarIndex]<0.5) {
						data.parameters.infoFile()<<"check track "<<trackCounter<<std::endl;
						data.parameters.infoFile()<<"Error: Inactive lifted edge on an active path"<<std::endl;
						data.parameters.infoFile().flush();
						std::cout<<"check track "<<trackCounter<<std::endl;
						throw std::runtime_error(std::string("Inactive lifted edge on an active path"));
					}
					else if(!isOnPath[vertex2]&&solution0[leVarIndex]>0.5){
						data.parameters.infoFile()<<"check track "<<trackCounter<<std::endl;
						data.parameters.infoFile()<<"Active lifted edge outside an active path"<<std::endl;
						data.parameters.infoFile().flush();
						std::cout<<"check track "<<trackCounter<<std::endl;
						throw std::runtime_error(std::string("Active lifted edge outside an active path"));
					}


				}


				size_t newVertex=0;
				bool nextFound=false;
				for (int j = 0; j < graph.numberOfEdgesToVertex(vertex);
						++j) {
					e = graph.edgeToVertex(vertex, j);
					eVarIndex = data.getEdgeVarIndex(e);

					if (solution0[eVarIndex] > 0.5) {
						newVertex = graph.vertexToVertex(vertex, j);
						path.push_front(newVertex);
						isOnPath[newVertex]=1;
						isOnAnyPath[newVertex]=1;
						activeEdgeCounter++;
						if(nextFound){
							data.parameters.infoFile()<<"check track "<<trackCounter<<std::endl;
							data.parameters.infoFile()<<"Multiple active incoming edges to vertex "<<vertex<<std::endl;
							data.parameters.infoFile().flush();
							std::cout<<"check track "<<trackCounter<<std::endl;
							std::cout<<"Multiple active incoming edges to vertex "<<vertex<<std::endl;
						}
						nextFound=true;


						//break;
					}
				}
				vertex=newVertex;



			}

			//std::cout<<std::endl;
			while(!path.empty()){
				isOnPath[path.front()]=0;
				path.pop_front();
			}
		}
	}

	for (int i = 0; i < graph.numberOfVertices(); ++i) {
		if(i==s||i==t) continue;
		if((solution0[i]>0.5&&!isOnAnyPath[i])||(solution0[i]<0.5&&isOnAnyPath[i])){
			data.parameters.infoFile()<<"check track "<<trackCounter<<std::endl;
			data.parameters.infoFile()<<"vertex "<<i<<" labeled "<<solution0[i]<<" is on Path "<<isOnAnyPath[i]<<std::endl;
			data.parameters.infoFile().flush();
			std::cout<<"check track "<<trackCounter<<std::endl;
			std::cout<<"vertex "<<i<<" labeled "<<solution0[i]<<" is on Path "<<isOnAnyPath[i]<<std::endl;
		}
	}
	size_t activeEdgeCounter2=0;
	for (int i = 0; i < graph.numberOfEdges(); ++i) {
		size_t index=data.getEdgeVarIndex(i);
		size_t v=data.getVertexVarIndex(graph.vertexOfEdge(i,0));
		size_t w=data.getVertexVarIndex(graph.vertexOfEdge(i,1));

		if(v==s||w==t) continue;
		if(solution0[index]>0.5){
			activeEdgeCounter2++;
			//if(v==s||w==t) continue;
			if(solution0[v]<0.5||solution0[w]<0.5){
				data.parameters.infoFile()<<"check track "<<trackCounter<<std::endl;
				data.parameters.infoFile()<<"edge "<<v<<","<<w<<" has label 1, but "<<solution0[v]<<","<<solution0[w]<<" are its vertices"<<std::endl;
				data.parameters.infoFile().flush();
				std::cout<<"check track "<<trackCounter<<std::endl;
				std::cout<<"edge "<<v<<","<<w<<" has label 1, but "<<solution0[v]<<","<<solution0[w]<<" are its vertices"<<std::endl;
			}
		}
		//		if(solution0[index]==0&&(solution0[v]==1||solution0[w]==1)){
		//			std::cout<<"edge "<<v<<","<<w<<" has label , but "<<solution0[v]<<","<<solution0[w]<<" are its vertices"<<std::endl;
		//		}
	}
	for (int i = 0; i < liftedGraph.numberOfEdges(); ++i) {
		size_t index=data.getLiftedEdgeVarIndex(i);
		size_t v=data.getVertexVarIndex(liftedGraph.vertexOfEdge(i,0));
		size_t w=data.getVertexVarIndex(liftedGraph.vertexOfEdge(i,1));
		if(v==s||w==t) continue;
		if(solution0[index]>0.5&&(solution0[v]<0.5||solution0[w]<0.5)){
			data.parameters.infoFile()<<"check track "<<trackCounter<<std::endl;
			data.parameters.infoFile()<<"lifted edge "<<v<<","<<w<<" has label 1, but "<<solution0[v]<<","<<solution0[w]<<" are its vertices"<<std::endl;
			data.parameters.infoFile().flush();
			std::cout<<"check track "<<trackCounter<<std::endl;
			std::cout<<"lifted edge "<<v<<","<<w<<" has label 1, but "<<solution0[v]<<","<<solution0[w]<<" are its vertices"<<std::endl;
		}
		//		if(solution0[index]==0&&(solution0[v]==1||solution0[w]==1)){
		//			std::cout<<"lifted edge "<<v<<","<<w<<" has label 0, but "<<solution0[v]<<","<<solution0[w]<<" are its vertices"<<std::endl;
		//		}
	}

	if(outputs){
		std::cout<<"Solution checked. Tracks: "<<trackCounter<<std::endl;
		data.parameters.infoFile()<<"Solution checked. Tracks: "<<trackCounter<<std::endl;
		data.parameters.infoFile().flush();
	}

}

template<class T=size_t>
void ilp_check_solution(Data<>& data,ilp::Gurobi& ilp0,bool outputs=true)
{
	std::vector<double> solution0(data.getCosts().size(),0);
	for (int i = 0; i < solution0.size(); ++i) {
		solution0[i]=ilp0.label(i);
		//if(ilp0.label(i)>0.5) std::cout<<ilp0.label(i)<<std::endl;
	}
	ilp_check_solution(data,solution0,outputs);
}



template<class T=size_t>
std::vector<double> ilp_solve(Data<>& data,bool outputs=true)
{
	ilp::Gurobi ilp;
	using Callback = disjointPaths::ilp::MyCallback;
	std::vector<double> initSolution;

	ilp.setVerbosity(true);
	if(outputs){
		if(data.getTask()=='T'){
			ilp.setLogFile(data.parameters.getOutputFileName() +"-tr-gurobi.log");
		}
		else{
			ilp.setLogFile(data.parameters.getOutputFileName() +"-gurobi.log");
		}
	}
	else{
		ilp.setOutputFlag(0);
		//ilp.setLogFile(data.parameters.getOutputFileName() +"-gurobi-temp.log");
	}
	ilp.setRelativeGap(data.parameters.getGurobiRelativeGap());


	ilp.setAbsoluteGap(0.0);

	ilp.addVariables(data.getCosts().size(),data.getCosts().data());

	//Callback callback(ilp, data,data.problemGraph.useExperimental_);
	Callback callback(ilp, data,outputs);
	//Callback callback(ilp, data,true);

	ilp.setCallback(callback);

	if(outputs){
		std::cout<<"callback set"<<std::endl;
		data.parameters.infoFile()<<"callback set"<<std::endl;
		data.parameters.infoFile().flush();
	}
	ilp_add_constraints<>(ilp,data);

	if(outputs){
		std::cout<<"constraints added"<<std::endl;
		data.parameters.infoFile()<<"constraints added"<<std::endl;
		data.parameters.infoFile().flush();
	}

	if(data.getInitSolution().size()==data.getCosts().size()){
		std::cout<<"Init solution from data structure"<<std::endl;
		data.parameters.infoFile()<<"Init solution from data structure"<<std::endl;
		data.parameters.infoFile().flush();
		initSolution=data.getInitSolution();
		ilp_check_solution(data,initSolution,true);
	}
	else{
		initSolution=ilp_init_solution<T>(data,outputs);
	}


	if(outputs){
		std::cout<<"set start"<<std::endl;
		data.parameters.infoFile()<<"set start"<<std::endl;
		data.parameters.infoFile().flush();
	}
	ilp.setStart(initSolution.begin());

	if(outputs){
		std::cout<<"opt start"<<std::endl;
		data.parameters.infoFile()<<"opt start"<<std::endl;
		data.parameters.infoFile().flush();
	}
	ilp.optimize();
	ilp_check_solution(data,ilp,outputs);

	std::vector<double> labels(data.getCosts().size());
	for (int i = 0; i < data.getCosts().size(); ++i) {
		double l=ilp.label(i);
		labels[i]=l;

	}
	return labels;

}




template<class T=size_t>
std::unordered_map<size_t,std::vector<size_t>> doPathOptimization(Data<>& data,std::vector<std::vector<size_t>>& allPaths)
{
	std::unordered_map<size_t,std::vector<size_t>> breaks;

	std::cout<<"all paths size "<<allPaths.size()<<std::endl;
	data.parameters.infoFile()<<"all paths size "<<allPaths.size()<<std::endl;
	for (int i = 0; i < allPaths.size(); ++i) {
		//std::cout<<"path "<<i<<" size :"<<allPaths[i].size()<<std::endl;
		//data.parameters.infoFile()<<"path "<<i<<" size :"<<allPaths[i].size()<<std::endl;
		if(allPaths[i].size()>1){
			data.onePathGraph(allPaths[i],false);
			//TODO init with all ones
			std::vector<double> labels=ilp_solve(data,false);
			//std::vector<std::vector<size_t>> fullPaths=data.pathsFromSolution(labels,false);
			std::vector<size_t> breaksInPath=data.cutsInOnePathSolution(labels);
			if(breaksInPath.size()>0){
				//TODO output path id and number of breaks
				std::cout<<"path "<<i<<" size :"<<allPaths[i].size()<<std::endl;
				data.parameters.infoFile()<<"path "<<i<<" size :"<<allPaths[i].size()<<std::endl;
				std::cout<<"breaks size "<<breaksInPath.size()<<std::endl;
				data.parameters.infoFile()<<"breaks size "<<breaksInPath.size()<<std::endl;
				breaks[i]=breaksInPath;
				//				std::cout<<"full paths size "<<fullPaths.size()<<std::endl;
				//				data.parameters.infoFile()<<"full paths size "<<fullPaths.size()<<std::endl;
			}
		}
	}
	std::cout<<"paths optimized, breaks size "<<breaks.size()<<std::endl;
	data.parameters.infoFile()<<"paths optimized, breaks size "<<breaks.size()<<std::endl;
	return breaks;

}


template<class T=size_t>
void solver_flow_only(ConfigDisjoint<>& parameters)
{
	std::cout<<"Using simple flow solver"<<std::endl;
	DisjointStructure<> DS(parameters);
	Data<> data(DS);
	std::vector<double> solution=ilp_init_solution(data,true);
	std::vector<std::vector<size_t>> paths=data.pathsFromSolution(solution,false);
	data.outputSolution(paths,false);
	std::cout<<"Simple flow solved"<<std::endl;
}


template<class T=size_t>
void solver_ilp(ConfigDisjoint<>& parameters)
{

	char delim=',';
	levinkov::Timer timer;
	timer.start();

	if(parameters.getMaxTimeLifted()==0){
		solver_flow_only(parameters);
	}
	else{

		CompleteStructure<> cs=CompleteStructure<>(parameters,delim);
		std::cout<<"vg max time "<<cs.vg.getMaxTime()<<std::endl;
		std::cout<<"vg max vertex "<<cs.vg.getMaxVertex()<<std::endl;
		std::cout<<"vg last frame size "<<cs.vg.getGroupVertices(cs.vg.getMaxTime()).size()<<std::endl;


		size_t maxTime=cs.maxTime;  //maxTime withou t
		size_t intervalSize=parameters.getSmallIntervals();
		size_t numberOfIntervals=maxTime/intervalSize;
		std::cout<<"cs max vertex: "<<cs.vg.getMaxVertex()<<", frame: "<<cs.vg.getGroupIndex(cs.vg.getMaxVertex())<<std::endl;
		parameters.infoFile()<<"cs max vertex: "<<cs.vg.getMaxVertex()<<", frame: "<<cs.vg.getGroupIndex(cs.vg.getMaxVertex())<<std::endl;
		//parameters.output("cs max vertex: "+cs.vg.getMaxVertex()+", frame: "+cs.vg.getGroupIndex(cs.vg.getMaxVertex())+"\n");

		if(maxTime%intervalSize!=0){
			numberOfIntervals++;
		}

		//bool halfIntervals=parameters.isOverlappingIntervals();
		bool halfIntervals=false;
		if(halfIntervals) numberOfIntervals=2*numberOfIntervals-1;

		std::vector<std::vector<size_t>> allPaths;
		bool loadIntervals=(parameters.getTaskType()=='T');
		if(loadIntervals){
			std::string fileName=parameters.getIntervalFile();
			allPaths=readLines<size_t>(fileName,' ');
			std::cout<<"loaded "<<allPaths.size()<<" paths"<<std::endl;
			parameters.infoFile()<<"loaded "<<allPaths.size()<<" paths"<<std::endl;

			//			std::vector<size_t> path={1328,1350,1362,1383,1478,1519,1529,1535,1547,1555,1566,1592,1600,1617,1631,1638,1737,1750,1801,1831,1842};
			//			allPaths.push_back(path);

		}
		else{
			int overlap=parameters.getOverlappingIntervals();
			for (int i = 0; i < numberOfIntervals; ++i) {
				levinkov::Timer intervalTimer;
				intervalTimer.start();
				size_t minT=i*intervalSize+1;
				size_t maxT=minT+intervalSize;
				int minTOverlap=minT-overlap;
				int maxTOverlap=maxT+overlap;
				if(halfIntervals){
					minT=(i*intervalSize)/2+1;
					maxT=minT+intervalSize;
				}


				maxT=std::min(maxT,maxTime+1);

				maxTOverlap=std::min(maxTOverlap,int(maxTime+1));
				minTOverlap=std::max(1,minTOverlap);

				std::cout<<"maxT "<<maxT<<std::endl;
				std::cout<<"maxTOverlap "<<maxTOverlap<<std::endl;
				std::cout<<"maxTime "<<maxTime<<std::endl;

				//				std::cout<<"interval "<<minT<<","<<maxT<<std::endl;
				//				parameters.infoFile()<<"interval "<<minT<<","<<maxT<<std::endl;
				//				parameters.infoFile().flush();
				//parameters.output("interval "+minT+","+maxT+"\n");
				//DisjointStructure<> DS=DisjointStructure<>(parameters,delim,&cs,minT,maxT);
				DisjointStructure<> DS=DisjointStructure<>(parameters,delim,&cs,minTOverlap,maxTOverlap);
				if(DS.getTerminalNode()==1){
					std::cout<<"Empty interval "<<std::endl;
					parameters.infoFile()<<"Empty interval "<<std::endl;
					 //empty interval
				}
				else if(DS.getGraph().numberOfEdges()<=(DS.getGraph().numberOfVertices()-2)*2){
					//s,t edges only
					std::cout<<"One layer interval "<<std::endl;
					parameters.infoFile()<<"One layer interval "<<std::endl;
					std::vector<std::vector<size_t>> newPaths(DS.maxV-DS.minV+1);
					for (int j = DS.minV; j <=DS.maxV; ++j) {
						std::vector<size_t> oneVertexPath(1);
						oneVertexPath[0]=j;
						newPaths[j-DS.minV]=oneVertexPath;
					}
					allPaths.insert(allPaths.end(),newPaths.begin(),newPaths.end());

				}
				else{


					Data<> data(DS);

					std::cout<<"min vertex in interval: "<<DS.minV<<std::endl;
					parameters.infoFile()<<"min vertex in interval: "<<DS.minV<<std::endl;
					std::cout<<"max vertex in interval: "<<data.getGraph().numberOfVertices()-3+DS.minV<<std::endl;
					parameters.infoFile()<<"max vertex in interval: "<<data.getGraph().numberOfVertices()-3+DS.minV<<std::endl;
					parameters.infoFile().flush();



					if(i==numberOfIntervals-1){
						std::cout<<"max vertex in last interval: "<<data.getGraph().numberOfVertices()-3+DS.minV<<std::endl;
						parameters.infoFile()<<"max vertex in last interval: "<<data.getGraph().numberOfVertices()-3+DS.minV<<std::endl;
						parameters.infoFile().flush();
					}


					std::vector<double> labels=ilp_solve(data);

					std::vector<std::vector<size_t>> newPaths=data.pathsFromSolution(labels,false,DS.minV);
//					std::cout<<"new paths in solver"<<std::endl;
//					for(auto& path:newPaths){
//						for(size_t v:path){
//							std::cout<<v<<", ";
//						}
//						std::cout<<std::endl;
//					}
					if(overlap>0){
						newPaths=extractInnerPaths(cs.vg,newPaths,minT,maxT-1);
					}

					if(halfIntervals){
						//TODO extract inner parts
						size_t minTForExtract;
						size_t maxTForExtract;
						if(i==0){
							minTForExtract=0;
						}
						else{
							minTForExtract=(i*intervalSize)/2+intervalSize/4+1;

						}
						if(i==numberOfIntervals-1){
							maxTForExtract=maxT-1;
						}
						else{
							maxTForExtract=((i+1)*intervalSize)/2+intervalSize/4;
						}
						newPaths=extractInnerPaths(cs.vg,newPaths,minTForExtract,maxTForExtract);
					}

					allPaths.insert(allPaths.end(),newPaths.begin(),newPaths.end());
					parameters.infoFile().flush();
				}

				intervalTimer.stop();
				timer.stop();
				double seconds=timer.get_elapsed_seconds();
				timer.start();
				std::cout<<"interval complete"<<std::endl;
				std::cout<<"interval solved in time "<<intervalTimer.get_elapsed_seconds()<<std::endl;
				std::cout<<"time from start "<<seconds<<std::endl;
				parameters.infoFile()<<"interval complete"<<std::endl;
				parameters.infoFile()<<"interval solved in time "<<intervalTimer.get_elapsed_seconds()<<std::endl;
				parameters.infoFile()<<"time from start "<<seconds<<std::endl;
				parameters.infoFile().flush();


			}
		}

		levinkov::Timer trackletTimer;
		trackletTimer.start();

		std::cout<<"all paths solved "<<std::endl;
		parameters.infoFile()<<"all paths solved "<<std::endl;
		parameters.infoFile().flush();

		bool optimizePaths=parameters.isOptimizePaths();

		Data<> data(cs,parameters);
		if(!optimizePaths){
//			std::cout<<"all paths in solver"<<std::endl;
//			for(auto& path:allPaths){
//				for(size_t v:path){
//					std::cout<<v<<", ";
//				}
//				std::cout<<std::endl;
//			}
			data.prepareGraphFromIntervalsDense(allPaths);
			if(!loadIntervals){
				data.outputSolution(allPaths,true);
			}

		}
		else{

			std::unordered_map<size_t,std::vector<size_t>> breaks=doPathOptimization(data,allPaths);
			data.graphFromOptimizedPaths(allPaths,breaks,false);
		}




		std::cout<<"call solver on all paths"<<std::endl;
		parameters.infoFile()<<"call solver on all paths"<<std::endl;
		parameters.infoFile().flush();

		std::vector<double> labels=ilp_solve(data);

		std::cout<<"obtain paths from solution"<<std::endl;
		parameters.infoFile()<<"obtain paths from solution"<<std::endl;
		parameters.infoFile().flush();

		std::vector<std::vector<size_t>> paths=data.pathsFromSolution(labels,true);
		data.evaluate(paths);

		bool isNewGraph;
		if(!optimizePaths){
			std::cout<<"call graph from interval dense for check"<<std::endl;
			parameters.infoFile()<<"call graph from interval dense for check"<<std::endl;
			parameters.infoFile().flush();
			isNewGraph=data.prepareGraphFromIntervalsDense(paths,true);
		}
		else{
			std::cout<<"call optimize paths for check"<<std::endl;
			parameters.infoFile()<<"call optimize paths for check"<<std::endl;
			parameters.infoFile().flush();
			std::unordered_map<size_t,std::vector<size_t>> breaks=doPathOptimization(data,paths);
			isNewGraph=data.graphFromOptimizedPaths(paths,breaks,true);
		}
		if(isNewGraph){
			std::cout<<"new graph created"<<std::endl;
			parameters.infoFile()<<"new graph created"<<std::endl;
			parameters.infoFile().flush();
		}
		else{
			std::cout<<"new graph not created"<<std::endl;
			parameters.infoFile()<<"new graph not created"<<std::endl;
			parameters.infoFile().flush();
		}
		while(isNewGraph){
			std::cout<<"final graph not optimal, recomputing..."<<std::endl;
			parameters.infoFile()<<"final graph not optimal, recomputing..."<<std::endl;
			parameters.infoFile().flush();

			labels=ilp_solve(data);
			paths=data.pathsFromSolution(labels,true);
			data.evaluate(paths);
			if(!optimizePaths){
				isNewGraph=data.prepareGraphFromIntervalsDense(paths,true);
			}
			else{
				std::unordered_map<size_t,std::vector<size_t>> breaks=doPathOptimization(data,paths);
				isNewGraph=data.graphFromOptimizedPaths(paths,breaks,true);
			}
		}

		data.outputSolution(paths);

		timer.stop();
		trackletTimer.stop();


		std::cout<<"Solver complete"<<std::endl;
		std::cout<<"Last phase time "<<trackletTimer.get_elapsed_seconds()<<std::endl;
		std::cout<<"Time from start "<<timer.get_elapsed_seconds()<<std::endl;
		parameters.infoFile()<<"Solver complete"<<std::endl;
		parameters.infoFile()<<"Last phase time "<<trackletTimer.get_elapsed_seconds()<<std::endl;
		parameters.infoFile()<<"Time from start "<<timer.get_elapsed_seconds()<<std::endl;
		parameters.infoFile().flush();
	}

}


//template<class ILP, class CALLBACK = ilp::Callback<ILP>>
template<class T=size_t>
void solver_ilp(DisjointStructure<>& ds)
{


	levinkov::Timer timer;
	timer.start();

	Data<> data(ds);

	std::vector<double> labels=ilp_solve(data);

	std::vector<std::vector<size_t>> paths=data.pathsFromSolution(labels,false);

	//If one more iteration would be done with cut into tracklets even if time breaks are empty, solution can be improved
	if(data.parameters.isDenseTracklets()){
		bool isNewGraph=data.prepareGraphFromIntervalsDense(paths,true);
		while(isNewGraph){
			std::cout<<"final graph not optimal, recomputing..."<<std::endl;
			ds.parameters.infoFile()<<"final graph not optimal, recomputing..."<<std::endl;
			ds.parameters.infoFile().flush();

			labels=ilp_solve(data);
			paths=data.pathsFromSolution(labels,true);
			isNewGraph=data.prepareGraphFromIntervalsDense(paths,true);
		}
	}

	data.outputSolution(paths);

	timer.stop();

	std::cout<<"solver complete"<<std::endl;
	std::cout<<"time "<<timer.get_elapsed_seconds()<<std::endl;
	ds.parameters.infoFile()<<"solver complete"<<std::endl<<"time "<<timer.get_elapsed_seconds()<<std::endl;
	ds.parameters.infoFile().flush();

}




//TODO export this as a debug procedure
//if(i==0){
//				std::vector<size_t> veticesToCheck={292,381};
//				for (int j = 0; j < veticesToCheck.size(); ++j) {
//					size_t vToCheck=veticesToCheck[j];
//					if(vToCheck>=DS.minV&&vToCheck<DS.minV+DS.getGraph().numberOfVertices()){
//						size_t v=vToCheck;
//						vToCheck-=DS.minV;
//						std::cout<<"base graph: "<<std::endl;
//						for (int k = 0; k < data.getGraph().numberOfEdgesFromVertex(vToCheck); ++k) {
//							size_t wToCheck=data.getGraph().vertexFromVertex(vToCheck,k);
//							size_t edgeToCheck=data.getGraph().edgeFromVertex(vToCheck,k);
//							double score=data.getCosts()[data.getEdgeVarIndex(edgeToCheck)];
//							size_t w=wToCheck+DS.minV;
//							std::cout<<v<<","<<w<<": "<<score<<std::endl;
//						}
//						for (int k = 0; k < data.getGraph().numberOfEdgesToVertex(vToCheck); ++k) {
//							size_t wToCheck=data.getGraph().vertexToVertex(vToCheck,k);
//							size_t edgeToCheck=data.getGraph().edgeToVertex(vToCheck,k);
//							double score=data.getCosts()[data.getEdgeVarIndex(edgeToCheck)];
//							size_t w=wToCheck+DS.minV;
//							std::cout<<w<<","<<v<<": "<<score<<std::endl;
//						}
//						std::cout<<"lifted graph: "<<std::endl;
//						for (int k = 0; k < data.getGraphLifted().numberOfEdgesFromVertex(vToCheck); ++k) {
//							size_t wToCheck=data.getGraphLifted().vertexFromVertex(vToCheck,k);
//							size_t edgeToCheck=data.getGraphLifted().edgeFromVertex(vToCheck,k);
//							double score=data.getCosts()[data.getLiftedEdgeVarIndex(edgeToCheck)];
//							size_t w=wToCheck+DS.minV;
//							std::cout<<v<<","<<w<<": "<<score<<std::endl;
//						}
//						for (int k = 0; k < data.getGraphLifted().numberOfEdgesToVertex(vToCheck); ++k) {
//							size_t wToCheck=data.getGraphLifted().vertexToVertex(vToCheck,k);
//							size_t edgeToCheck=data.getGraphLifted().edgeToVertex(vToCheck,k);
//							double score=data.getCosts()[data.getLiftedEdgeVarIndex(edgeToCheck)];
//							size_t w=wToCheck+DS.minV;
//							std::cout<<w<<","<<v<<": "<<score<<std::endl;
//						}
//					}
//				}
//}





} // namespace disjoint-paths



#endif /* INCLUDE_DISJOINT_PATHS_ILP_SOLVER_DISJOINT_ILP_HXX_ */
