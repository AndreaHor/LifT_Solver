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
#include "disjoint-paths/disjointPathsData.hxx"
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
        data.parameters.getControlOutput()<<"init solution"<<std::endl;
        data.parameters.writeControlOutput();
    }
    std::vector<double> solution0(data.getCosts().size(),0);
    for (int i = 0; i < solution0.size(); ++i) {
        solution0[i]=ilp0.label(i);
        //if(ilp0.label(i)>0.5) std::cout<<ilp0.label(i)<<std::endl;
    }

    double costDifference=0;


    if(outputs){
        data.parameters.getControlOutput()<<"fix init solution"<<std::endl;
        data.parameters.writeControlOutput();
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

                            data.parameters.getControlOutput()<<"Init: Two active incoming edges to vertex "<<std::endl;
                            data.parameters.writeControlOutput();
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
        data.parameters.getControlOutput()<<"Init solution: inactive on path: "<<inactiveOnPath<<", active outside path: "<<activeOutsidePath<<std::endl;
        data.parameters.getControlOutput()<<"cost difference "<<costDifference<<std::endl;
        data.parameters.writeControlOutput();
    }
    return solution0;



}



template<class T=size_t>
void ilp_check_solution(Data<>& data,std::vector<double>& solution0,bool outputs=true)
{

    if(outputs){
        data.parameters.getControlOutput()<<"fix first solution"<<std::endl;
        data.parameters.writeControlOutput();
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
                        data.parameters.getControlOutput()<<"check track "<<trackCounter<<std::endl;
                        data.parameters.getControlOutput()<<"Error: Inactive lifted edge on an active path"<<std::endl;
                        data.parameters.writeControlOutput();
                        throw std::runtime_error(std::string("Inactive lifted edge on an active path"));
                    }
                    else if(!isOnPath[vertex2]&&solution0[leVarIndex]>0.5){
                        data.parameters.getControlOutput()<<"check track "<<trackCounter<<std::endl;
                        data.parameters.getControlOutput()<<"Active lifted edge outside an active path"<<std::endl;
                        data.parameters.writeControlOutput();
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
                            data.parameters.getControlOutput()<<"check track "<<trackCounter<<std::endl;
                            data.parameters.getControlOutput()<<"Multiple active incoming edges to vertex "<<vertex<<std::endl;
                            data.parameters.writeControlOutput();

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
            data.parameters.getControlOutput()<<"check track "<<trackCounter<<std::endl;
            data.parameters.getControlOutput()<<"vertex "<<i<<" labeled "<<solution0[i]<<" is on Path "<<isOnAnyPath[i]<<std::endl;
            data.parameters.writeControlOutput();

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
                data.parameters.getControlOutput()<<"check track "<<trackCounter<<std::endl;
                data.parameters.getControlOutput()<<"edge "<<v<<","<<w<<" has label 1, but "<<solution0[v]<<","<<solution0[w]<<" are its vertices"<<std::endl;
                data.parameters.writeControlOutput();

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
            data.parameters.getControlOutput()<<"check track "<<trackCounter<<std::endl;
            data.parameters.getControlOutput()<<"lifted edge "<<v<<","<<w<<" has label 1, but "<<solution0[v]<<","<<solution0[w]<<" are its vertices"<<std::endl;
            data.parameters.writeControlOutput();

        }
        //		if(solution0[index]==0&&(solution0[v]==1||solution0[w]==1)){
        //			std::cout<<"lifted edge "<<v<<","<<w<<" has label 0, but "<<solution0[v]<<","<<solution0[w]<<" are its vertices"<<std::endl;
        //		}
    }

    if(outputs){
        data.parameters.getControlOutput()<<"Solution checked. Tracks: "<<trackCounter<<std::endl;
        data.parameters.writeControlOutput();
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

        data.parameters.getControlOutput()<<"callback set"<<std::endl;
        data.parameters.writeControlOutput();
    }
    ilp_add_constraints<>(ilp,data);

    if(outputs){
        data.parameters.getControlOutput()<<"constraints added"<<std::endl;
        data.parameters.writeControlOutput();
    }

    if(data.getInitSolution().size()==data.getCosts().size()){
        data.parameters.getControlOutput()<<"Init solution from data structure"<<std::endl;
        data.parameters.writeControlOutput();
        initSolution=data.getInitSolution();
        ilp_check_solution(data,initSolution,true);
    }
    else{
        initSolution=ilp_init_solution<T>(data,outputs);
    }


    if(outputs){

        data.parameters.getControlOutput()<<"set start"<<std::endl;
        data.parameters.writeControlOutput();
    }
    ilp.setStart(initSolution.begin());

    if(outputs){

        data.parameters.getControlOutput()<<"opt start"<<std::endl;
        data.parameters.writeControlOutput();
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

    data.parameters.getControlOutput()<<"all paths size "<<allPaths.size()<<std::endl;
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
                data.parameters.getControlOutput()<<"path "<<i<<" size :"<<allPaths[i].size()<<std::endl;
                data.parameters.getControlOutput()<<"breaks size "<<breaksInPath.size()<<std::endl;
                breaks[i]=breaksInPath;
                //				std::cout<<"full paths size "<<fullPaths.size()<<std::endl;
                //				data.parameters.infoFile()<<"full paths size "<<fullPaths.size()<<std::endl;
            }
        }
    }

    data.parameters.getControlOutput()<<"paths optimized, breaks size "<<breaks.size()<<std::endl;
    data.parameters.writeControlOutput();
    return breaks;

}


template<class T=size_t>
std::vector<std::vector<size_t>> solver_flow_only(DisjointParams<>& parameters)
{
   parameters.getControlOutput()<<"Using simple flow solver"<<std::endl;
    DisjointStructure<> DS(parameters);
    Data<> data(DS);
    std::vector<double> solution=ilp_init_solution(data,true);
    std::vector<std::vector<size_t>> paths=data.pathsFromSolution(solution,false);
    //data.outputSolution(paths,false);
    parameters.getControlOutput()<<"Simple flow solved"<<std::endl;
    parameters.writeControlOutput();
    return paths;
}

template<class T=size_t>
std::vector<std::vector<size_t>> solver_ilp_tracklets(DisjointParams<>& parameters,CompleteStructure<>& cs,  std::vector<std::vector<size_t>>& allPaths,bool loadIntervals)
{
    levinkov::Timer trackletTimer;
    trackletTimer.start();

    parameters.getControlOutput()<<"all paths solved "<<std::endl;
    parameters.writeControlOutput();

    bool optimizePaths=parameters.isOptimizePaths();

    Data<> data(cs,parameters);
    if(!optimizePaths){

        data.prepareGraphFromIntervalsDense(allPaths);
        if(!loadIntervals&&parameters.isControlOutputFilesSet()){
            data.outputSolution(allPaths,true);
        }

    }
    else{

        std::unordered_map<size_t,std::vector<size_t>> breaks=doPathOptimization(data,allPaths);
        data.graphFromOptimizedPaths(allPaths,breaks,false);
    }




    parameters.getControlOutput()<<"call solver on all paths"<<std::endl;
    parameters.writeControlOutput();

    std::vector<double> labels=ilp_solve(data,parameters.isControlOutputFilesSet());

    parameters.getControlOutput()<<"obtain paths from solution"<<std::endl;
    parameters.writeControlOutput();

    std::vector<std::vector<size_t>> paths=data.pathsFromSolution(labels,true);
    data.evaluate(paths);

    bool isNewGraph;
    if(!optimizePaths){
        parameters.getControlOutput()<<"call graph from interval dense for check"<<std::endl;
        parameters.writeControlOutput();
        isNewGraph=data.prepareGraphFromIntervalsDense(paths,true);
    }
    else{

        parameters.getControlOutput()<<"call optimize paths for check"<<std::endl;
        parameters.writeControlOutput();
        std::unordered_map<size_t,std::vector<size_t>> breaks=doPathOptimization(data,paths);
        isNewGraph=data.graphFromOptimizedPaths(paths,breaks,true);
    }
    if(isNewGraph){

        parameters.getControlOutput()<<"new graph created"<<std::endl;
        parameters.writeControlOutput();
    }
    else{

        parameters.getControlOutput()<<"new graph not created"<<std::endl;
        parameters.writeControlOutput();
    }
    while(isNewGraph){

        parameters.getControlOutput()<<"final graph not optimal, recomputing..."<<std::endl;
        parameters.writeControlOutput();

        labels=ilp_solve(data,parameters.isControlOutputFilesSet());
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

   // data.outputSolution(paths);

    trackletTimer.stop();

    parameters.getControlOutput()<<"Solver complete"<<std::endl;
    parameters.getControlOutput()<<"Last phase time "<<trackletTimer.get_elapsed_seconds()<<std::endl;
    parameters.writeControlOutput();

    return paths;


}

template<class T=size_t>
std::vector<std::vector<size_t>> solver_ilp_no_intervals(DisjointParams<>& parameters,CompleteStructure<>& cs){
    levinkov::Timer timer;
    timer.start();

    size_t minT=1;
    size_t maxT=cs.maxTime+1;  //plus one because DS solver excludes maxT
    char delim=',';
    DisjointStructure<> DS=DisjointStructure<>(parameters,delim,&cs,minT,maxT);

    Data<> data0(DS);


    parameters.getControlOutput()<<"min vertex: "<<DS.getMinV()<<std::endl;
    parameters.getControlOutput()<<"max vertex: "<<data0.getGraph().numberOfVertices()-3+DS.getMinV()<<std::endl;
    parameters.writeControlOutput();

    std::vector<double> labels=ilp_solve(data0,parameters.isControlOutputFilesSet());

    std::vector<std::vector<size_t>> newPaths=data0.pathsFromSolution(labels,false);

    if(parameters.isDenseTracklets()){
        std::vector<std::vector<size_t>> paths=solver_ilp_tracklets(parameters,cs,newPaths,false);
        timer.stop();

        parameters.getControlOutput()<<"Time from start "<<timer.get_elapsed_seconds()<<std::endl;
        parameters.writeControlOutput();
        return paths;
    }
    else{
     //   data0.outputSolution(newPaths);
        timer.stop();

        parameters.getControlOutput()<<"Time from start "<<timer.get_elapsed_seconds()<<std::endl;
        parameters.writeControlOutput();
        return newPaths;
    }


}




template<class T=size_t>
std::vector<std::vector<size_t>> solver_ilp_intervals(DisjointParams<>& parameters,CompleteStructure<>& cs)
{

    char delim=',';
    levinkov::Timer timer;
    timer.start();

    std::vector<std::vector<size_t>> paths;

    if(parameters.getMaxTimeLifted()==0){
        paths=solver_flow_only(parameters);
    }
    else{

        //CompleteStructure<> cs=CompleteStructure<>(parameters,delim);
//        std::cout<<"vg max time "<<cs.getVertexGroups().getMaxTime()<<std::endl;
//        std::cout<<"vg max vertex "<<cs.getVertexGroups().getMaxVertex()<<std::endl;
//        std::cout<<"vg last frame size "<<cs.getVertexGroups().getGroupVertices(cs.getVertexGroups().getMaxTime()).size()<<std::endl;

        parameters.getControlOutput()<<"vg max time "<<cs.getVertexGroups().getMaxTime()<<std::endl;
        parameters.getControlOutput()<<"vg max vertex "<<cs.getVertexGroups().getMaxVertex()<<std::endl;
        parameters.getControlOutput()<<"vg last frame size "<<cs.getVertexGroups().getGroupVertices(cs.getVertexGroups().getMaxTime()).size()<<std::endl;
        parameters.writeControlOutput();

        size_t maxTime=cs.maxTime;  //maxTime withou t
        size_t intervalSize=parameters.getSmallIntervals();
        size_t numberOfIntervals=maxTime/intervalSize;


//        std::cout<<"cs max vertex: "<<cs.getVertexGroups().getMaxVertex()<<", frame: "<<cs.getVertexGroups().getGroupIndex(cs.getVertexGroups().getMaxVertex())<<std::endl;
//        parameters.infoFile()<<"cs max vertex: "<<cs.getVertexGroups().getMaxVertex()<<", frame: "<<cs.getVertexGroups().getGroupIndex(cs.getVertexGroups().getMaxVertex())<<std::endl;


        parameters.getControlOutput()<<"cs max vertex: "<<cs.getVertexGroups().getMaxVertex()<<", frame: "<<cs.getVertexGroups().getGroupIndex(cs.getVertexGroups().getMaxVertex())<<std::endl;
        parameters.writeControlOutput();


        //parameters.output("cs max vertex: "+cs.vg.getMaxVertex()+", frame: "+cs.vg.getGroupIndex(cs.vg.getMaxVertex())+"\n");

        if(maxTime%intervalSize!=0){
            numberOfIntervals++;
        }


        //bool halfIntervals=false;
        std::vector<std::vector<size_t>> allPaths;
        bool loadIntervals=(parameters.getTaskType()=='T');
        if(loadIntervals){
            std::string fileName=parameters.getIntervalFile();
            allPaths=readLines<size_t>(fileName,' ');
            parameters.getControlOutput()<<"loaded "<<allPaths.size()<<" paths"<<std::endl;

            parameters.writeControlOutput();

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

                maxT=std::min(maxT,maxTime+1);

                maxTOverlap=std::min(maxTOverlap,int(maxTime+1));
                minTOverlap=std::max(1,minTOverlap);

                parameters.getControlOutput()<<"maxT "<<maxT<<std::endl;
                parameters.getControlOutput()<<"maxTOverlap "<<maxTOverlap<<std::endl;
                parameters.getControlOutput()<<"maxTime "<<maxTime<<std::endl;
                parameters.writeControlOutput();


                DisjointStructure<> DS=DisjointStructure<>(parameters,delim,&cs,minTOverlap,maxTOverlap);  //last time used maxTOverlap-1
                if(DS.getTerminalNode()==1){
                    levinkov::Timer timer;
                    parameters.getControlOutput()<<"Empty interval "<<std::endl;
                    parameters.writeControlOutput();
                }
                else if(DS.getGraph().numberOfEdges()<=(DS.getGraph().numberOfVertices()-2)*2){
                    //s,t edges only
                    parameters.getControlOutput()<<"One layer interval "<<std::endl;
                    parameters.writeControlOutput();
                    std::vector<std::vector<size_t>> newPaths(DS.getMaxV()-DS.getMinV()+1);
                    for (int j = DS.getMinV(); j <=DS.getMaxV(); ++j) {
                        std::vector<size_t> oneVertexPath(1);
                        oneVertexPath[0]=j;
                        newPaths[j-DS.getMinV()]=oneVertexPath;
                    }
                    allPaths.insert(allPaths.end(),newPaths.begin(),newPaths.end());

                }
                else{


                    Data<> data(DS);

                    parameters.getControlOutput()<<"min vertex in interval: "<<DS.getMinV()<<std::endl;
                    parameters.getControlOutput()<<"max vertex in interval: "<<data.getGraph().numberOfVertices()-3+DS.getMinV()<<std::endl;
                    parameters.writeControlOutput();



                    if(i==numberOfIntervals-1){
                        parameters.getControlOutput()<<"max vertex in last interval: "<<data.getGraph().numberOfVertices()-3+DS.getMinV()<<std::endl;
                        parameters.writeControlOutput();
                    }


                    std::vector<double> labels=ilp_solve(data,parameters.isControlOutputFilesSet());

                    std::vector<std::vector<size_t>> newPaths=data.pathsFromSolution(labels,false,DS.getMinV());

                    if(overlap>0){
                        newPaths=cs.getVertexGroups().extractInnerPaths(newPaths,minT,maxT-1);
                    }


                    allPaths.insert(allPaths.end(),newPaths.begin(),newPaths.end());
                    parameters.writeControlOutput();
                }

                intervalTimer.stop();
                timer.stop();
                double seconds=timer.get_elapsed_seconds();
                timer.start();
                parameters.getControlOutput()<<"interval complete"<<std::endl;
                parameters.getControlOutput()<<"interval solved in time "<<intervalTimer.get_elapsed_seconds()<<std::endl;
                parameters.getControlOutput()<<"time from start "<<seconds<<std::endl;
                parameters.writeControlOutput();


            }
        }
        paths=solver_ilp_tracklets<>(parameters,cs,allPaths,loadIntervals);

        timer.stop();
        parameters.getControlOutput()<<"Time from start "<<timer.get_elapsed_seconds()<<std::endl;
        parameters.writeControlOutput();

    }
    return paths;
}




template<class T=size_t>
std::vector<std::vector<size_t>> solver_ilp(DisjointParams<>& parameters,CompleteStructure<>& cs)
{
    if(parameters.getSmallIntervals()>0){
        std::vector<std::vector<size_t>> paths=solver_ilp_intervals(parameters,cs);
        return paths;
    }
    else{
        std::vector<std::vector<size_t>> paths=solver_ilp_no_intervals(parameters,cs);
        return paths;
    }

}







} // namespace disjointPaths



#endif /* INCLUDE_DISJOINT_PATHS_ILP_SOLVER_DISJOINT_ILP_HXX_ */
