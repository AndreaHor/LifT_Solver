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
#include "disjoint-paths/parametersParser.hxx"


struct Parameters {
	std::string configFile;

};


Parameters parseCommandLine(int argc, char** argv){
try
{
    Parameters parameters;

      TCLAP::CmdLine tclap("track", ' ', "1.0");


    TCLAP::ValueArg<std::string> argConfigFile("c", "congig-file", "congig-file", true, parameters.configFile, "congig-file", tclap);


    tclap.parse(argc, argv);
    parameters.configFile = argConfigFile.getValue();


return parameters;
}
catch(TCLAP::ArgException& e)
{
    throw std::runtime_error(e.error());
}

}


int main(int argc, char** argv)
try
{

	auto parameters = parseCommandLine(argc, argv);
    disjointPaths::ParametersParser parametersParser;
    parametersParser.initFromFile(parameters.configFile);

    disjointPaths::DisjointParams<size_t> disjointParameters(parametersParser.getParsedStrings());
    if(disjointParameters.isParametersSet()){


            disjointPaths::CompleteStructure<> cs(disjointParameters);
            std::vector<std::vector<size_t>> paths=disjointPaths::solver_ilp<size_t>(disjointParameters,cs);
            disjointPaths::writeOutputToFile(paths,disjointParameters.getOutputFileName() + "-all-paths-FINAL.txt");


	}
	else{
		std::cout<<"Problem with parameter settings. No computation done."<<std::endl;
	}




    return 0;
}
catch (const std::runtime_error& error)
{
    std::cerr << "error: " << error.what() << std::endl;
    return 1;
}
