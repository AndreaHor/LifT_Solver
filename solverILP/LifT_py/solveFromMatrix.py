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

#Initializes structure for holding solver parameters. It expects the path to the solver parameter file.
params=ldpPy.DisjointParams("../data/exampleSolverILP/params_sequence.ini")

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

#Calling se solver on the given problem. 
ldpPy.solve_ilp(params,completeGraphStructure)
	
