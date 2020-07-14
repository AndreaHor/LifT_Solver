#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 11:37:09 2020

@author: fuksova
"""


import disjointPathsPy as ldpPy

params=ldpPy.DisjointParams("../data/exampleSolverILP/params_sequence.ini")

timeFrames=ldpPy.TimeFramesToVertices()
timeFrames.init_from_file("../data/exampleSolverILP/problemDesc_frames",params)

completeGraphStructure=ldpPy.GraphStructure(timeFrames)
completeGraphStructure.add_edges_from_file("../data/exampleSolverILP/problemDesc",params)
ldpPy.solve_ilp(params,completeGraphStructure)
	
