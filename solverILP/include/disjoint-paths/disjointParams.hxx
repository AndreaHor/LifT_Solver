/*
 * inputConfix.hxx
 *
 *  Created on: Jun 27, 2019
 *      Author: fuksova
 */

#ifndef INCLUDE_DISJOINT_PATHS_INPUTCONFIG_HXX_
#define INCLUDE_DISJOINT_PATHS_INPUTCONFIG_HXX_

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <limits>
#include"disjoint-paths/parametersParser.hxx"

namespace disjointPaths {


//template<class T=size_t>
//std::vector<std::string> split(const
//                               std::string& inputString, char delim) {
//    size_t occurence = 0;
//    size_t newOccurence = 0;
//    std::vector<std::string> strings;
//    while (newOccurence < inputString.size()) {
//        newOccurence = std::min(inputString.find_first_of(delim, occurence),
//                                inputString.size());

//        std::string newString(inputString, occurence, newOccurence - occurence);
//        strings.push_back(newString);
//        newOccurence = newOccurence + 1;
//        occurence = newOccurence;
//    }

//    return strings;
//}


template<class T = size_t>
class DisjointParams {
public:

    ~DisjointParams(){
        getControlOutput()<<"config disjoint destructor"<<std::endl;
        writeControlOutput();
        if(pInfoFile!=nullptr){
            (*pInfoFile).close();
            delete pInfoFile;
        }
    }


   // DisjointParams(const std::string& fileName);
    // DisjointParams(ParametersParser& paramsConstructor);
     DisjointParams(std::map<std::string,std::string>& parsedParams);
    //DisjointParamsDisjointParamsDisjointParams(const std::stringstream& stream);

    template<class STR>
    std::map<std::string,std::string> parseStream(STR& data);

    const std::string& getGraphFileName() const {
        return graphFileName;
    }

    size_t getMaxTimeBase() const {
        return maxTimeBase;
    }

    size_t getMaxTimeFrame() const {
        return maxTimeFrame;
    }

    size_t getMaxTimeLifted() const {
        return maxTimeLifted;
    }

    const std::string& getOutputFileName() const {
        return outputFileName;
    }

    const std::string& getTimeFileName() const {
        return timeFileName;
    }

    bool isAutomaticLifted() const {
        return automaticLifted;
    }

    double getBaseUpperThreshold() const {
        return baseUpperThreshold;
    }

    size_t getDenseTimeLifted() const {
        return denseTimeLifted;
    }

    size_t getKnnK() const {
        return knnK;
    }

    size_t getKnnTimeGap() const {
        return knnTimeGap;
    }

    size_t getLongerIntervalLifted() const {
        return longerIntervalLifted;
    }

    double getNegativeThresholdLifted() const {
        return negativeThresholdLifted;
    }

    double getPositiveThresholdLifted() const {
        return positiveThresholdLifted;
    }

    bool isRepulsive() const {
        return repulsive;
    }

    bool isRestrictFrames() const {
        return restrictFrames;
    }

    double getInputCost() const {
        return inputCost;
    }

    double getOutputCost() const {
        return outputCost;
    }

    bool isSparsify() const {
        return sparsify;
    }

    size_t getSmallIntervals() const {
        return smallIntervals;
    }

    double getGurobiRelativeGap() const {
        return gurobiRelativeGap;
    }

    bool isRequireImproving() const {
        return requireImproving;
    }

    double getGurobiRelGapTrackletG() const {
        return gurobiRelGapTrackletG;
    }

    bool isDebugOutputFiles() const {
        return debugOutputFiles;
    }

    bool isNewLifted() const {
        return newLifted;
    }

    size_t getTrackletSize() const {
        return trackletSize;
    }

    bool isDenseTracklets() const {
        return denseTracklets;
    }

    void setMaxTimeLifted(size_t maxTimeLifted) {
        this->maxTimeLifted = maxTimeLifted;
    }

    int getOverlappingIntervals() const {
        return overlappingIntervals;
    }

    size_t getMaxTimeGapComplete() const {
        return maxTimeGapComplete;
    }

    //	void output(std::string toOutput){
    //		(*fileForGeneralOutputs)<<toOutput;
    //	}



    std::stringstream& getControlOutput(){  //TODO make private
        return controlOutput;
    }

    bool isAllBaseTracklet() const {
        return allBaseTracklet;
    }

    bool isOptimizePaths() const {
        return optimizePaths;
    }

    size_t getRepulsiveTimeGap() const {
        return repulsiveTimeGap;
    }

    bool isCompleteGapForTracklet() const {
        return completeGapForTracklet;
    }

    size_t getTrTimeGapBase() const {
        return trTimeGapBase;
    }

    size_t getTrTimeGapLifted() const {
        return trTimeGapLifted;
    }

    const std::string& getIntervalFile() const {
        return intervalFile;
    }

    char getTaskType() const {
        return taskType;
    }

    bool isParametersSet() const {
        return parametersSet;
    }

    //ConfigDisjoint<>& operator=(const ConfigDisjoint<>&);


    bool isControlOutputFilesSet() const{
        return controlOutputFiles;
    }

    void writeControlOutput(){
        if(stdOutput){
            std::cout<<controlOutput.str();

        }

        if(controlOutputFiles){
            infoFile()<<controlOutput.str();
            infoFile().flush();
        }
        controlOutput=std::stringstream();
    }

private:
    DisjointParams<T>(const DisjointParams<T>&);
   // void initParameters(std::map<std::string,std::string>& parameters);
    DisjointParams();

    std::ofstream& infoFile(){
        std::ofstream & infoF=*pInfoFile;
        return infoF;
    }
    //std::ofstream * fileForGeneralOutputs;
    std::ofstream* pInfoFile;
    std::string graphFileName;
    std::string outputFileName;
    std::string timeFileName;
    std::string paramsFileName;
    bool debugOutputFiles;
    bool controlOutputFiles;
    bool stdOutput;

    std::stringstream controlOutput;

    double inputCost;
    double outputCost;

    bool sparsify;
    bool restrictFrames;
    size_t maxTimeFrame;
    size_t smallIntervals;


    size_t maxTimeBase;
    double baseUpperThreshold;
    size_t knnTimeGap;
    size_t knnK;
    bool requireImproving;

    bool automaticLifted;
    double negativeThresholdLifted;
    double positiveThresholdLifted;
    size_t maxTimeLifted;
    size_t denseTimeLifted;
    size_t longerIntervalLifted;
    size_t repulsiveTimeGap;

    size_t trTimeGapLifted;
    size_t trTimeGapBase;

    size_t trackletSize;
    bool denseTracklets;
    int overlappingIntervals;
    bool optimizePaths;
    size_t maxTimeGapComplete;

    bool allBaseTracklet;
    bool newLifted;
    double gurobiRelativeGap;
    double gurobiRelGapTrackletG;

    bool completeGapForTracklet;
    char taskType;
    std::string intervalFile;

    bool parametersSet;


    bool repulsive;

};


// template<class T>
// inline DisjointParams<T>::DisjointParams(const std::stringstream& stream){
// 	std::map<std::string,std::string> params=parseStream(stream);
// 	initParameters(params);
// }


//template<class T>
//inline DisjointParams<T>::DisjointParams(const std::string& fileName){
//    parametersSet=false;
//    pInfoFile=nullptr;
//    std::ifstream data;
//    //data.exceptions( std::ifstream::failbit | std::ifstream::badbit );
//    try{
//        data.open(fileName);
//        if (!data){
//            throw std::system_error(errno, std::system_category(), "failed to open "+fileName);
//        }
//        std::map<std::string,std::string> params=parseStream(data);
//        initParameters(params);

//    }

//    catch (std::system_error& er) {
//        std::clog << er.what() << " (" << er.code() << ")" << std::endl;
//    }

//}

//template<class T>
//inline DisjointParams<T>::DisjointParams(ParametersParser& paramsConstructor){
//    DisjointParams(paramsConstructor.getParsedStrings());
//}

template<class T>
inline DisjointParams<T>::DisjointParams(std::map<std::string,std::string>& parameters){
    pInfoFile=nullptr;
    //std::map<std::string,std::string>& parameters=paramsConstructor.getParsedStrings();
    if(parameters.count("CONTROL_OUTPUT_FILES")>0){
        controlOutputFiles=std::stoi(parameters["CONTROL_OUTPUT_FILES"]);
    }
    else{
        controlOutputFiles=true;
    }
    if(parameters.count("CONTROL_STD_OUTPUT")>0){
        stdOutput=std::stoi(parameters["CONTROL_OUTPUT_FILES"]);
    }
    else{
        stdOutput=true;
    }

    controlOutput<<"control output files "<<controlOutputFiles<<std::endl;




    if(parameters.count("INPUT_GRAPH")>0){
        graphFileName=parameters["INPUT_GRAPH"];
    }
    else{
        graphFileName="";
        //throw std::runtime_error("File with the input graph was not specified in the config file");
    }
    controlOutput<<"input graph "<<graphFileName<<std::endl;

    std::string pathToInputFile=graphFileName.substr(0,graphFileName.find_last_of('/')+1);

    if(parameters.count("OUTPUT_PATH")==0){
        if(graphFileName.empty()){
            parameters["OUTPUT_PATH"]="./";
        }
        else{
            parameters["OUTPUT_PATH"]=pathToInputFile;
        }

    }

    //TODO fix output prefix, it is surrounded by spaces!
    if(parameters.count("OUTPUT_PREFIX")==0){
        parameters["OUTPUT_PREFIX"]="output";
    }


    outputFileName=parameters["OUTPUT_PATH"]+parameters["OUTPUT_PREFIX"];
    paramsFileName=outputFileName+"-params.txt";
    std::ofstream paramsFile;
    if(controlOutputFiles){
         paramsFile.open(paramsFileName.data(),std::ofstream::out);
    }


    if(controlOutputFiles){
        controlOutput<<"params file "<<paramsFileName<<std::endl;
        paramsFile<<"control output files "<<controlOutputFiles<<std::endl;
        paramsFile<<"input graph "<<graphFileName<<std::endl;
        paramsFile<<"output files "<<outputFileName<<std::endl;
        paramsFile<<"params file "<<paramsFileName<<std::endl;
    }
    std::string infoFileName=outputFileName+"-info.txt";

    if(controlOutputFiles){
        pInfoFile=new std::ofstream(infoFileName.data(),std::ofstream::out);  //TODO only if controlOutputFiles  is set
    }
    writeControlOutput();

    //fileForGeneralOutputs=&paramsFile;

    //TODO use params File name to output all params displayed here with cout

    controlOutput<<"output files "<<outputFileName<<std::endl;

    if(parameters.count("INPUT_FRAMES")>0){
        timeFileName=parameters["INPUT_FRAMES"];
    }
    else{
        timeFileName=graphFileName+"_frames";
    }
    controlOutput<<"input frames "<<timeFileName<<std::endl;
     if(controlOutputFiles) paramsFile<<"input frames "<<timeFileName<<std::endl;


    if(parameters.count("DEBUG_OUTPUT_FILES")>0){
        debugOutputFiles=std::stoi(parameters["DEBUG_OUTPUT_FILES"]);
    }
    else{
        debugOutputFiles=0;
    }
    controlOutput<<"debug output files "<<debugOutputFiles<<std::endl;
    if(controlOutputFiles) paramsFile<<"debug output files "<<debugOutputFiles<<std::endl;



    if(parameters.count("SPARSIFY")>0){
        sparsify=std::stoi(parameters["SPARSIFY"]);
    }
    else{
        sparsify=1;
    }
    controlOutput<<"sparsify "<<sparsify<<std::endl;
    if(controlOutputFiles) paramsFile<<"sparsify "<<sparsify<<std::endl;

    //TODO Put some of the following parameters under if(sparsify){...
    //TODO check for wrong input?
    //TODO cout or write the params unto output file or both

    if(parameters.count("RESTRICT_FRAMES")>0){
        restrictFrames=std::stoi(parameters["RESTRICT_FRAMES"]);
    }
    else{
        restrictFrames=1;
    }
    controlOutput<<"restrict frames "<<restrictFrames<<std::endl;
    if(controlOutputFiles) paramsFile<<"restrict frames "<<restrictFrames<<std::endl;

    if(parameters.count("MAX_TIMEGAP")>0){
        maxTimeFrame=std::stoul(parameters["MAX_TIMEGAP"]);
        if(maxTimeFrame==0){maxTimeFrame=std::numeric_limits<size_t>::max();}
    }
    else{
        maxTimeFrame=std::numeric_limits<size_t>::max();
    }
    std::cout<<"max time gap "<<maxTimeFrame<<std::endl;
    if(controlOutputFiles) paramsFile<<"max time gap "<<maxTimeFrame<<std::endl;

    if(parameters.count("INPUT_COST")>0){
        inputCost=std::stod(parameters["INPUT_COST"]);
    }
    else{
        inputCost=0;
    }
    controlOutput<<"input cost "<<inputCost<<std::endl;
    if(controlOutputFiles) paramsFile<<"input cost "<<inputCost<<std::endl;

    if(parameters.count("OUTPUT_COST")>0){
        outputCost=std::stod(parameters["OUTPUT_COST"]);
    }
    else{
        outputCost=0;
    }
    controlOutput<<"output cost "<<outputCost<<std::endl;
    if(controlOutputFiles) paramsFile<<"output cost "<<outputCost<<std::endl;

    writeControlOutput();

    if(parameters.count("MAX_TIMEGAP_BASE")>0){
        maxTimeBase=std::stoul(parameters["MAX_TIMEGAP_BASE"]);
    }
    else{
        maxTimeBase=60;
    }
    controlOutput<<"max time gap base "<<maxTimeBase<<std::endl;
    if(controlOutputFiles) paramsFile<<"max time gap base "<<maxTimeBase<<std::endl;

    if(parameters.count("KNN_GAP")>0){
        knnTimeGap=std::stoul(parameters["KNN_GAP"]);
    }
    else{
        knnTimeGap=3;
    }
    controlOutput<<"KNN time gap "<<knnTimeGap<<std::endl;
    if(controlOutputFiles) paramsFile<<"KNN time gap "<<knnTimeGap<<std::endl;

    if(parameters.count("KNN_K")>0){
        knnK=std::stoul(parameters["KNN_K"]);
    }
    else{
        knnK=3;
    }
    controlOutput<<"KNN k "<<knnK<<std::endl;
    if(controlOutputFiles) paramsFile<<"KNN k "<<knnK<<std::endl;

    if(parameters.count("BASE_THRESHOLD")>0){
        baseUpperThreshold=std::stod(parameters["BASE_THRESHOLD"]);
    }
    else{
        baseUpperThreshold=0;
    }
    controlOutput<<"base upper threshold "<<baseUpperThreshold<<std::endl;
    if(controlOutputFiles) paramsFile<<"base upper threshold "<<baseUpperThreshold<<std::endl;

    writeControlOutput();
    //	if(parameters.count("REQUIRE_IMPROVING")>0){
    //		requireImproving=std::stoi(parameters["REQUIRE_IMPROVING"]);
    //	}
    //	else{
    requireImproving=0;
    //}
    controlOutput<<"require improving "<<requireImproving<<std::endl;
    if(controlOutputFiles) paramsFile<<"require improving "<<requireImproving<<std::endl;

    //	if(parameters.count("AUTOMATIC_LIFTED")>0){
    //		automaticLifted=std::stoi(parameters["AUTOMATIC_LIFTED"]);
    //	}
    //	else{
    automaticLifted=true;
    //}
    controlOutput<<"automatic lifted "<<automaticLifted<<std::endl;
    if(controlOutputFiles) paramsFile<<"automatic lifted "<<automaticLifted<<std::endl;

    if(parameters.count("MAX_TIMEGAP_LIFTED")>0){
        maxTimeLifted=std::stoul(parameters["MAX_TIMEGAP_LIFTED"]);
    }
    else{
        maxTimeLifted=60;
    }
    controlOutput<<"max time gap lifted "<<maxTimeLifted<<std::endl;
    if(controlOutputFiles) paramsFile<<"max time gap lifted "<<maxTimeLifted<<std::endl;

    if(parameters.count("DENSE_TIMEGAP_LIFTED")>0){
        denseTimeLifted=std::stoul(parameters["DENSE_TIMEGAP_LIFTED"]);
    }
    else{
        denseTimeLifted=12;
    }
    controlOutput<<"dense time gap lifted "<<denseTimeLifted<<std::endl;
    if(controlOutputFiles) paramsFile<<"dense time gap lifted "<<denseTimeLifted<<std::endl;

    if(parameters.count("NEGATIVE_THRESHOLD_LIFTED")>0){
        negativeThresholdLifted=std::stod(parameters["NEGATIVE_THRESHOLD_LIFTED"]);
    }
    else{
        negativeThresholdLifted=-1;
    }
    controlOutput<<"negative threshold lifted "<<negativeThresholdLifted<<std::endl;
    if(controlOutputFiles) paramsFile<<"negative threshold lifted "<<negativeThresholdLifted<<std::endl;

    if(parameters.count("POSITIVE_THRESHOLD_LIFTED")>0){
        positiveThresholdLifted=std::stod(parameters["POSITIVE_THRESHOLD_LIFTED"]);
    }
    else{
        positiveThresholdLifted=1;
    }
    controlOutput<<"positive threshold lifted "<<positiveThresholdLifted<<std::endl;
    if(controlOutputFiles) paramsFile<<"positive threshold lifted "<<positiveThresholdLifted<<std::endl;

    if(parameters.count("LONGER_LIFTED_INTERVAL")>0){
        longerIntervalLifted=std::stoul(parameters["LONGER_LIFTED_INTERVAL"]);
    }
    else{
        longerIntervalLifted=4;
    }
    controlOutput<<"longer interval lifted "<<longerIntervalLifted<<std::endl;
    if(controlOutputFiles) paramsFile<<"longer interval lifted "<<longerIntervalLifted<<std::endl;


    //	if(parameters.count("REPULSIVE_TIMEGAP")>0){
    //		repulsiveTimeGap=std::stoul(parameters["REPULSIVE_TIMEGAP"]);
    //	}
    //	else{
    repulsiveTimeGap=maxTimeLifted;
    //}
    //std::cout<<"repulsive time gap "<<repulsiveTimeGap<<std::endl;
    //paramsFile<<"repulsive time gap "<<repulsiveTimeGap<<std::endl;


    //	if(parameters.count("NEW_LIFTED")>0){
    //		newLifted=std::stoi(parameters["NEW_LIFTED"]);
    //	}
    //	else{
    newLifted=false;
    //	}
    //	std::cout<<"new lifted sparse. strategy "<<newLifted<<std::endl;
    //	paramsFile<<"new lifted sparse. strategy "<<newLifted<<std::endl;




    if(parameters.count("SMALL_INTERVALS")>0){
        smallIntervals=std::stoul(parameters["SMALL_INTERVALS"]);
    }
    else{
        smallIntervals=40;
    }
    controlOutput<<"small intervals "<<smallIntervals<<std::endl;
    if(controlOutputFiles) paramsFile<<"small intervals "<<smallIntervals<<std::endl;


    writeControlOutput();


    if(parameters.count("TRACKLET_SIZE")>0){
        trackletSize=std::stoul(parameters["TRACKLET_SIZE"]);
        if(smallIntervals>0){
            trackletSize=std::min(trackletSize,smallIntervals);
        }
    }
    else{
        trackletSize=20;
        trackletSize=std::min(trackletSize,smallIntervals);
    }
    controlOutput<<"tracklet size "<<trackletSize<<std::endl;
    if(controlOutputFiles) paramsFile<<"tracklet size "<<trackletSize<<std::endl;

    if(parameters.count("ALL_BASE_TRACKLET")>0){
        allBaseTracklet=std::stoi(parameters["ALL_BASE_TRACKLET"]);

    }
    else{
        allBaseTracklet=true;
    }
    controlOutput<<"all base edges in tracklet graph "<<allBaseTracklet<<std::endl;
    if(controlOutputFiles) paramsFile<<"all base edges in tracklet graph "<<allBaseTracklet<<std::endl;



    if(parameters.count("DENSE_TRACKLETS")>0){
        denseTracklets=std::stoi(parameters["DENSE_TRACKLETS"]);
    }
    else{
        denseTracklets=true;
    }
    controlOutput<<"dense tracklets "<<denseTracklets<<std::endl;
    if(controlOutputFiles) paramsFile<<"dense tracklets "<<denseTracklets<<std::endl;


    if(parameters.count("OPTIMIZE_PATHS")>0){
        optimizePaths=std::stoi(parameters["OPTIMIZE_PATHS"]);
    }
    else{
        optimizePaths=false;
    }
    controlOutput<<"optimize paths "<<optimizePaths<<std::endl;
    if(controlOutputFiles) paramsFile<<"optimize paths "<<optimizePaths<<std::endl;

    if(parameters.count("MAX_TIMEGAP_COMPLETE")>0){
        maxTimeGapComplete=std::stoul(parameters["MAX_TIMEGAP_COMPLETE"]);
        maxTimeGapComplete=std::max(maxTimeGapComplete,maxTimeBase);
        maxTimeGapComplete=std::max(maxTimeGapComplete,maxTimeLifted);
    }
    else{
        maxTimeGapComplete=std::max(maxTimeLifted,maxTimeBase);
    }
    controlOutput<<"max time gap complete "<<maxTimeGapComplete<<std::endl;
    if(controlOutputFiles) paramsFile<<"max time gap complete "<<maxTimeGapComplete<<std::endl;



    if(parameters.count("COMPLETE_GAP_IN_TRACKLET")>0){
        completeGapForTracklet=std::stoi(parameters["COMPLETE_GAP_IN_TRACKLET"]);
    }
    else{
        completeGapForTracklet=false;
    }
    controlOutput<<"use complete time gap in tracklet graph "<<completeGapForTracklet<<std::endl;
    if(controlOutputFiles) paramsFile<<"use complete time gap in tracklet graph "<<completeGapForTracklet<<std::endl;



    //	if(parameters.count("OVERLAP_INTERVALS")>0){
    //		overlappingIntervals=std::stoi(parameters["OVERLAP_INTERVALS"]);
    //		if(smallIntervals==0) overlappingIntervals=0;
    //	}
    //	else{
    overlappingIntervals=0;
    //}
    //std::cout<<"overlapping intervals "<<overlappingIntervals<<std::endl;
    //paramsFile<<"overlapping intervals "<<overlappingIntervals<<std::endl;

    //
    //	if(parameters.count("REPULSIVE")>0){
    //		repulsive=std::stoi(parameters["REPULSIVE"]);
    //	}
    //	else{
    repulsive=0;
    //	}
    controlOutput<<"use repulsive "<<repulsive<<std::endl;
    if(controlOutputFiles) paramsFile<<"use repulsive "<<repulsive<<std::endl;

    if(parameters.count("GUROBI_REL_GAP")>0){
        gurobiRelativeGap=std::stod(parameters["GUROBI_REL_GAP"]);
    }
    else{
        gurobiRelativeGap=0;
    }
    controlOutput<<"gurobi relative gap "<<gurobiRelativeGap<<std::endl;
    if(controlOutputFiles) paramsFile<<"gurobi relative gap "<<gurobiRelativeGap<<std::endl;


    if(parameters.count("GUROBI_REL_GAP_TRACKLET")>0){
        gurobiRelGapTrackletG=std::stod(parameters["GUROBI_REL_GAP_TRACKLET"]);
    }
    else{
        gurobiRelGapTrackletG=0;
    }
    controlOutput<<"gurobi relative gap for tracklet graph "<<gurobiRelGapTrackletG<<std::endl;
    if(controlOutputFiles) paramsFile<<"gurobi relative gap for tracklet graph "<<gurobiRelGapTrackletG<<std::endl;


    if(parameters.count("TASK")>0){
        std::string task=parameters["TASK"];
        taskType=task[0];
        //std::cout<<"loaded letter"<<taskType<<std::endl;
        if(taskType!='T'&&taskType!='C'){
            taskType='C';
        }
    }
    else{
        taskType='C';
    }
    controlOutput<<"task type "<<taskType<<std::endl;
    if(controlOutputFiles) paramsFile<<"task type "<<taskType<<std::endl;

    if(taskType=='T'){
        if(parameters.count("INTERVAL_FILE")>0){
            intervalFile=parameters["INTERVAL_FILE"];
        }
        else{
            intervalFile=parameters["OUTPUT_PATH"]+parameters["OUTPUT_PREFIX"]+"-all-paths-INTERVALS.txt";
        }
        controlOutput<<"Interval solution from "<<intervalFile<<std::endl;
        if(controlOutputFiles) paramsFile<<"Interval solution from "<<intervalFile<<std::endl;
    }
    else{
        intervalFile="";
    }


    trTimeGapBase=getMaxTimeGapComplete();
    trTimeGapLifted=getMaxTimeGapComplete();

    if(!isCompleteGapForTracklet()){
        trTimeGapBase=getMaxTimeBase();
        trTimeGapLifted=getMaxTimeLifted();
    }

    controlOutput<<"max time gap base for tracklet graph "<<trTimeGapBase<<std::endl;
    if(controlOutputFiles) paramsFile<<"max time gap base for tracklet graph "<<trTimeGapBase<<std::endl;


    controlOutput<<"max time gap lifted for tracklet graph "<<trTimeGapLifted<<std::endl;
    if(controlOutputFiles) paramsFile<<"max time gap lifted for tracklet graph "<<trTimeGapLifted<<std::endl;


    if(controlOutputFiles) paramsFile.close();
    parametersSet=true;

}




//template<class T>
//template<class STR>
//inline std::map<std::string,std::string> DisjointParams<T>::parseStream(STR& data){
//    char delim='=';
//    std::string line;
//    std::vector<std::string> strings;
//    std::map<std::string,std::string> parameters;



//    size_t lineCounter=0;
//    std::vector<size_t> currentGroup;
//    bool solverPart=false;
//    while (!solverPart&&std::getline(data, line) ) {
//        size_t commentPos=line.find_first_of('#');
//        if(commentPos<line.length()){
//            line.erase(commentPos);
//        }
//        if(!line.empty()&&line.find("[SOLVER]")!=line.npos){
//            solverPart=true;
//        }
//    }
//    if(!solverPart){
//        throw std::runtime_error("Config file does not contain \"[SOLVER]\".  ");
//    }
//    bool newSection=false;
//    while(!newSection&&std::getline(data, line)){
//        size_t commentPos=line.find_first_of('#');
//        if(commentPos<line.length()){
//            line.erase(commentPos);
//        }
//        if(line.find("[")==line.npos){
//            if(!line.empty()){ //TODO not split all delims, just the first occurence
//                strings=split<>(line,delim);
//                std::string whitespaces (" ");

//                size_t foundLast = strings[0].find_last_not_of(whitespaces);
//                size_t foundFirst=strings[0].find_first_not_of(whitespaces);
//                std::string key=strings[0].substr(foundFirst,foundLast-foundFirst+1);

//                foundLast = strings[1].find_last_not_of(whitespaces);
//                foundFirst=strings[1].find_first_not_of(whitespaces);
//                std::string value=strings[1].substr(foundFirst,foundLast-foundFirst+1);

//                parameters[key]=value;
//            }
//        }
//        else{
//            newSection=true;
//        }
//    }

//    for(auto itMap=parameters.begin();itMap!=parameters.end();itMap++){
//        std::cout<<itMap->first<<"-"<<itMap->second<<"-"<<std::endl;
//    }
//    return parameters;
//}










}//End namespace disjointPaths

#endif /* INCLUDE_DISJOINT_PATHS_INPUTCONFIG_HXX_ */
