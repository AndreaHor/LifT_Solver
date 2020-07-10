#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 11:37:09 2020

@author: fuksova
"""

#import numpy as np
#import disjointPathsPy as ldpPy

params=ldpPy.DisjointParams("/home/fuksova/codes/higher-order-disjoint-paths/params_test.ini")

timeFrames=ldpPy.TimeFramesToVertices()
timeFrames.init_from_file("/home/fuksova/codes/higher-order-disjoint-paths/data/error13_5_20/problemDesc_frames",params)

completeGraphStructure=ldpPy.GraphStructure(timeFrames)
completeGraphStructure.add_edges_from_file("/home/fuksova/codes/higher-order-disjoint-paths/data/error13_5_20/problemDesc",params)
ldpPy.solve_ilp(completeGraphStructure)

