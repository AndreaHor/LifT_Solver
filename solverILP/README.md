# LifT_Solver
Solver for lifted disjoint paths problem

This is a c++ solver for the lifted disjoint paths problem. Given an input graph structure with edge costs, it extracts the base graph and the lifted graph from it and solves the lifted disjoint paths problem on those graphs. Its output is the list of found tracks. Each track is a list of graph vertices.

The solver can either be used independently or within our tracker pipeline. The whole pipeline for the tracker is contained in the directory LifT_Solver/tracker. The tracker provides suitable input files for the c++ solver, runs the solver on them and processes its output.

A documentation for the c++ solver itself is contained in LifT_Solver/solverILP/documentation/solverDocumentation.pdf.
