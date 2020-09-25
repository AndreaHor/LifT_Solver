#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 11:37:09 2020

@author: fuksova
"""

#This is a script showing how to use LifT ILP solver from python. Its usage is similar to directly calling run-disjoint-paths 
#from command line. The solver parameters are read from a parameter file.
# Graph structure together with assignemt of vertices to time frames are initialized from the same two files that are used for
#direct calling of run-disjoint-paths. 
#The main difference is that the paths to the two files containing graph structure are not read from the parameter file.
#They have to be provided as in the script.

import disjointPathsPy as ldpPy

#Create a parser for parameters
paramsParser=ldpPy.ParametersParser()

#Parses parameters from file
paramsParser.init_from_file("../data/exampleSolverILP/params_sequence_py.ini")

#Initializes structure for holding solver parameters. It expects a string to string map as an input. ParametersParser.get_parsed_params() can be used for providing such map.
params=ldpPy.DisjointParams(paramsParser.get_parsed_params())

#Constructor of structure for holding the mapping between time frames and graph vertices
timeFrames=ldpPy.TimeFramesToVertices()

#Initalizing the structure from a file
timeFrames.init_from_file("../data/exampleSolverILP/problemDesc_frames",params)

#Initializing the graph structure from timeFrames. For now, no edges are present. 
completeGraphStructure=ldpPy.GraphStructure(timeFrames)

#Adding edges to graph structure from a file.
completeGraphStructure.add_edges_from_file("../data/exampleSolverILP/problemDesc",params)

#Calling the solver on the given problem. 
paths=ldpPy.solve_ilp(params,completeGraphStructure)

#The resulting paths are allways automatically saved into the file with suffix "-all-paths-FINAL" in the output directory. 
#You can optionally save the resulting paths into another file. 
ldpPy.write_output_to_file(paths,"../data/exampleSolverILP/my_python_output.txt")

for path in paths:
  for v in path:
    print(v, end =" ")
  print("")
