#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 11:37:09 2020

@author: fuksova
"""


import disjointPathsPy as ldpPy
import numpy as np

params=ldpPy.DisjointParams("../data/exampleSolverILP/params_sequence_py.ini")

timeFrames=ldpPy.TimeFramesToVertices()
vect= np.array([3, 2, 3])
timeFrames.init_from_vector(vect)

completeGraphStructure=ldpPy.GraphStructure(timeFrames)
matrix12=np.array([[-2,-1],[-1,-2],[-1,-1]])
matrix23=np.array([[-2,-1,-1],[-1,-2,-1]])
matrix13=np.array([[-2,-1,-1],[-1,-2,-1],[-1,-1,-2]])
completeGraphStructure.add_edges_from_array(1,2,matrix12)
completeGraphStructure.add_edges_from_array(2,3,matrix23)
completeGraphStructure.add_edges_from_array(1,3,matrix13)

ldpPy.solve_ilp(params,completeGraphStructure)
	
