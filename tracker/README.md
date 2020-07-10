# Lif_T Tracker 

## Setup

1. Make sure gurobi and the Lift solver are installed 

2. Adapt all paths in runtime_configuration/configuration.ini:
	1. `datasetsDir = LifT_Solver/data/tmp/` points to the root folder that contains different MOT datasets
	2. `featuresDir = LifT_Solver/tracker/features/data/tmp/` points to the folder containing the graph weights
	3. `solverDir = LifT_Solver/solverILP` points to the folder containing the lifted disjoint paths solver
	4. `outputDir = LifT_Solver/tracker/results/` points to the folder where the results are stored


