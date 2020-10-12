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
//#include "disjoint-paths/disjointStructure.hxx"
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



}
}




#endif /* INCLUDE_DISJOINT_PATHS_ILP_MYCALLBACK_HXX_ */
