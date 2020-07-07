/*
 * runDisjointPahts.cxx
 *
 *  Created on: Sep 10, 2018
 *      Author: fuksova
 */

#include <stdexcept>

//#include <andres/ilp/gurobi-callback.hxx>

#include <tclap/CmdLine.h>

#include "andres/graph/graph.hxx"
#include "disjoint-paths/disjointPathsData.hxx"
//#include "disjoint-paths/ilp/gurobi-callback-disjoint.hxx"
#include "disjoint-paths/ilp/solver-disjoint-ilp.hxx"
//#include "disjoint-paths/ilp/MyCallback.hxx"
#include "disjoint-paths/disjointParams.hxx"
#include "disjoint-paths/disjointPathsMethods.hxx"


struct Parameters {
	std::string configFile;
//    std::string graphFileName;
//    std::string timeFileName="";
//    bool useRepulsive=false;
//    bool useConnectivity=true;
//
//    bool useExperimental=true;
//    std::string outputName="";
//    size_t maxTimeBase=1;
//    size_t maxTimeLifted=0;
//    size_t maxTimeFrame=1000;
};


Parameters parseCommandLine(int argc, char** argv)
try
{
    Parameters parameters;

      TCLAP::CmdLine tclap("track", ' ', "1.0");
//    TCLAP::ValueArg<std::string> argGraphFileName("g", "graph-file", "graph information", true, parameters.graphFileName, "graph-file", tclap);
//    TCLAP::ValueArg<std::string> argTimeFileName("t", "time-file", "time information", false, parameters.timeFileName, "time-file", tclap);
//
//    TCLAP::ValueArg<std::string> argOutputName("o", "output-name", "output name", false, parameters.outputName, "output-name", tclap);
//
//    TCLAP::ValueArg<bool> argUseRepulsive("r", "use-repulsive", "use repulsive", false, parameters.useRepulsive, "use-repulsive", tclap);
//    TCLAP::ValueArg<bool> argUseExperimental("x", "use-experimental", "use experimental", false, parameters.useExperimental, "use-experimental", tclap);
//    TCLAP::ValueArg<bool> argUseConnectivity("c", "use-connectivity", "use connectivity", false, parameters.useConnectivity, "use-connectivity", tclap);
//    TCLAP::ValueArg<size_t> argMaxTimeBase("b", "max-time-base", "max-time-base", false, parameters.maxTimeBase, "max-time-base", tclap);
//    TCLAP::ValueArg<size_t> argMaxTimeLifted("l", "max-time-lifted", "max-time-lifted", false, parameters.maxTimeLifted, "max-time-lifted", tclap);
//    TCLAP::ValueArg<size_t> argMaxTimeFrame("m", "max-time-frame", "max-time-frame", false, parameters.maxTimeFrame, "max-time-frame", tclap);

    TCLAP::ValueArg<std::string> argConfigFile("c", "congig-file", "congig-file", true, parameters.configFile, "congig-file", tclap);


    tclap.parse(argc, argv);
    parameters.configFile = argConfigFile.getValue();

//    parameters.graphFileName = argGraphFileName.getValue();
//    parameters.useRepulsive=argUseRepulsive.getValue();
//    if(parameters.useRepulsive){
//    	parameters.useRepulsive=false;
//    	std::cout<<"Warning: Repulsive costs currently unsupported! Running version without repulsive costs."<<std::endl;
//    }
//    parameters.outputName=argOutputName.getValue();
//    parameters.useConnectivity=argUseConnectivity.getValue();
//    parameters.useExperimental=argUseExperimental.getValue();
//    parameters.timeFileName=argTimeFileName.getValue();
//    if(!parameters.timeFileName.empty()){
//    	parameters.maxTimeBase=argMaxTimeBase.getValue();
//    	parameters.maxTimeLifted=argMaxTimeLifted.getValue();
//    	parameters.maxTimeFrame=argMaxTimeFrame.getValue();
//    }
//
//    if(parameters.outputName.empty()){
//    	size_t indexMax=parameters.graphFileName.find_last_of(".");
//    	size_t indexMin=parameters.graphFileName.find_last_of("/");
//    	if((indexMin<parameters.graphFileName.npos&&indexMin<indexMax)||indexMin==std::string::npos){
//    		parameters.outputName=parameters.graphFileName.substr(0,indexMax);
//    	}
//    	else{
//    		parameters.outputName=parameters.graphFileName;
//    	}
//    }

    return parameters;
}
catch(TCLAP::ArgException& e)
{
    throw std::runtime_error(e.error());
}


int main(int argc, char** argv)
try
{

	//disjointPaths::ConfigDisjoint<> config("/home/fuksova/codes/higher-order-disjoint-paths/matlab/Tracking/runtime_configuration/paramsBigger.ini");


	auto parameters = parseCommandLine(argc, argv);
	disjointPaths::DisjointParams<size_t> configDisjoint(parameters.configFile);
	if(configDisjoint.isParametersSet()){
		if(configDisjoint.getMaxTimeLifted()==0){
			disjointPaths::solver_flow_only(configDisjoint);
		}
		else{

			if(configDisjoint.getSmallIntervals()==0){
				disjointPaths::DisjointStructure<> disjointP(configDisjoint);
				disjointPaths::solver_ilp<size_t>(disjointP);
			}
			else{
				disjointPaths::CompleteStructure<> cs(configDisjoint);
				cs.addEdgesFromFile();
				disjointPaths::solver_ilp<size_t>(configDisjoint,cs);
			}
		}
	}
	else{
		std::cout<<"Problem with parameter settings. No computation done."<<std::endl;
	}


//	auto parameters = parseCommandLine(argc, argv);
//	if(parameters.timeFileName.empty()){
//		disjointPaths::DisjointStructure<> disjointP(parameters.graphFileName,',',parameters.useExperimental);
//		disjointPaths::solver_ilp<disjointPaths::ilp::Gurobi>(disjointP,parameters.outputName,parameters.useRepulsive);
//	}
//	else{
//		disjointPaths::DisjointStructure<> disjointP(parameters.graphFileName,parameters.timeFileName,parameters.maxTimeBase,parameters.maxTimeLifted,parameters.maxTimeFrame,',');
//		disjointPaths::solver_ilp<disjointPaths::ilp::Gurobi>(disjointP,parameters.outputName,parameters.useRepulsive);
//	}

	// disjointPaths::DisjointInit<> disjointP("input.txt");


    return 0;
}
catch (const std::runtime_error& error)
{
    std::cerr << "error: " << error.what() << std::endl;
    return 1;
}
