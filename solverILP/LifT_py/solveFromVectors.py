#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 11:37:09 2020

@author: fuksova
"""

#This is a script showing how to use LifT ILP solver from python. Here, solver parameters are given in a python dictionary.
#Graph structure together with assignemt of vertices to time frames are initialized from python arrays.

import disjointPathsPy as ldpPy
import numpy as np


paramsMap={}
paramsMap["INPUT_COST"]="0"
paramsMap["OUTPUT_PREFIX"]="config_output_"
paramsMap["INPUT_COST"]="0"
paramsMap["OUTPUT_COST"]="0"
paramsMap["CONF_COST"]="0"
paramsMap["GREEDY"]="0"
paramsMap["GREEDY_SIZE"]="1"
paramsMap["LIFTED_SECONDS"]="2"
paramsMap["BASE_SECONDS"]="2"
paramsMap["KNN_GAP"]="3"
paramsMap["KNN_K"]="3"
paramsMap["BASE_THRESHOLD"]="-0.5"
paramsMap["DENSE_TIMEGAP_LIFTED"]="12"
paramsMap["NEGATIVE_THRESHOLD_LIFTED_OLD"]="-0.3"
paramsMap["NEGATIVE_THRESHOLD_LIFTED"]="-2.2"
paramsMap["POSITIVE_THRESHOLD_LIFTED_OLD"]="0.3"
paramsMap["POSITIVE_THRESHOLD_LIFTED"]="2.2"
paramsMap["LONGER_LIFTED_INTERVAL"]="4"
paramsMap["SMALL_INTERVALS"]="40"
paramsMap["TRACKLET_SIZE"]="20"
paramsMap["GUROBI_REL_GAP"]="0"
paramsMap["GUROBI_REL_GAP_TRACKLET"]="0"
paramsMap["DEBUG_OUTPUT_FILES"]="0"
paramsMap["CONTROL_OUTPUT_FILES"]="0"
paramsMap["ALL_BASE_TRACKLET"]="1"
paramsMap["OPTIMIZE_PATHS"]="0"
paramsMap["MAX_TIMEGAP_BASE"]="60"
paramsMap["MAX_TIMEGAP_LIFTED"]="60"
paramsMap["OUTPUT_PATH"]="../data/exampleSolverILP/"


#Initializes structure for holding solver parameters. It expects a string to string map (dictionary) as an input. ParametersParser.get_parsed_params() can be alternatively used for providing such map.
params=ldpPy.DisjointParams(paramsMap)

#Constructor of structure for holding the mapping between time frames and graph vertices
timeFrames=ldpPy.TimeFramesToVertices()

#This array contains information about how many detections are in each time frame. The array must contain non-negative integer values. 
#Numbering of graph vertices in time frames: Frame 1: 0,1,2; Frame 2: 3,4; Frame 3: 4,5,6
vect= np.array([3, 2, 3])

#Initalizing the structure from the given array.
timeFrames.init_from_vector(vect)

#Initializing the graph structure from timeFrames. For now, no edges are present. 
completeGraphStructure=ldpPy.GraphStructure(timeFrames)

#Edge vector. Each edge is represented by a pair of vertices. Vertices ar numbered from zero
edgeVector=np.array([[0,3],[0,4],[1,3],[1,4],[2,3],[2,4],[3,5],[3,6],[3,7],[4,5],[4,6],[4,7],[0,5],[0,6],[0,7],[1,5],[1,6],[1,7],[2,6],[2,7]])
costVector=np.array([-2.2,-1,-1,-2.2,-1,-1,-2.2,-1,-1,-1,-2.2,-1,  -2.2,-1,-1, -1,-2.2,-1, -1,-2.2])

#Adding edges to graph structure from vectors
completeGraphStructure.add_edges_from_vectors(edgeVector,costVector,params)

#Calling the solver on the given problem and get the resulting paths. 
paths=ldpPy.solve_ilp(params,completeGraphStructure)

#Saving the resulting paths into an output file. 
ldpPy.write_output_to_file(paths,"../data/exampleSolverILP/my_python_output_vect.txt")

for path in paths:
  for v in path:
    print(v, end =" ")
  print("")
