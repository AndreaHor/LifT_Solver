#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 11:37:09 2020

@author: fuksova
"""

#This is a script showing how to use LifT ILP solver from python. Here, solver parameters are read from
#a parameter file.
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


#Initializes structure for holding solver parameters. It expects either a ParameterParser or a string to string dictionary as an input
params=ldpPy.DisjointParams(paramsMap)

#Constructor of structure for holding the mapping between time frames and graph vertices
timeFrames=ldpPy.TimeFramesToVertices()

#This array contains information about how many detections are in each time frame. The array must contain non-negative integer values. 
vect= np.array([3, 2, 3])

#Initalizing the structure from the given array.
timeFrames.init_from_vector(vect)

#Initializing the graph structure from timeFrames. For now, no edges are present. 
completeGraphStructure=ldpPy.GraphStructure(timeFrames)


#This ininity variable will be used as an edge weight indicating that an edge is missing.
noEdge=float("inf")

#Three adjacency matrices that specify costs of edges between pairs of time frames. The number of rows of each matrix must correspond to the number of detections in the lower time frame. The number of matrix columns must correspond to the number of detections in the higher time frame.
matrix12=np.array([[-2.2,-1.0],[-1.0,-2.2],[-1.0,-1.0]])
matrix23=np.array([[-2.2,-1.0,-1.0],[-1.0,-2.2,-1.0]])
matrix13=np.array([[-2.2,-1.0,2.3],[-1.0,-2.2,-1.0],[noEdge,-1.0,-2.2]])

#Adding edges to graph structure from adjacency matrices.
#Arguments: index of the lower frame, index of the higher frame, adjacency matrix
completeGraphStructure.add_edges_from_array(1,2,matrix12)
completeGraphStructure.add_edges_from_array(2,3,matrix23)
completeGraphStructure.add_edges_from_array(1,3,matrix13)

#Calling the solver on the given problem and get the resulting paths. 
paths=ldpPy.solve_ilp(params,completeGraphStructure)

#The resulting paths are allways automatically saved into the file with suffix "-all-paths-FINAL" in the output directory. 
#You can optionally save the resulting paths into another file. 
ldpPy.write_output_to_file(paths,"../data/exampleSolverILP/my_python_output.txt")

for path in paths:
  for v in path:
    print(v, end =" ")
  print("")
