# Lif_T Tracker 

We provide a wrapper around the Lift solver that loads detections and graph weights, then calls the solver and stores the tracking results in the MOTChallenge format. It can also produce result videos. 

Not contained in this repository are the feature computation, pre- and post-processing.


## Setup

1. Make sure gurobi and the Lift solver are installed 


2. `runtime_configuration/configuration.ini`: Stores references to ressources, and parameters that can be used to achieve speed-ups for a parameter study.
	1. `datasetsDir = LifT_Solver/data/tmp/` points to the root folder that contains different MOT datasets
	2. `featuresDir = LifT_Solver/tracker/features/data/tmp/` points to the folder containing the graph weights
	3. `solverDir = LifT_Solver/solverILP` points to the folder containing the lifted disjoint paths solver
	4. `outputDir = LifT_Solver/tracker/results/` points to the folder where the results are stored
	5. `saveImg`: Boolean variable deciding wether to create result videos. Works only if FFMPEG is accessable
	6. `CALL_SOLVER`,`EVAL_SOLVER`, `CREATE_PROBLEM_FILE`: They can be used to speed-up the tracker if one sequence is computed with different parameters. If different solver parameters are used, `CREATE_PROBLEM_FILE` can be set to 0 after the first run. For a study on the post-processing, `CALL_SOLVER` and `CREATE_PROBLEM_FILE` can be set to 0. 
	7. `gurobiLic` points to the gurobi license file 

3. `runtime_configuration/params.ini`: Contains all relevant parameters for the solver as well as the tracker.
	1. For all parameters related to the solver, please look in the documentary of the solver.
	2. The parameters `LIFTED_SECONDS` and `BASE_SECONDS` control the maximal edge length in seconds.

4. Precomputed features are provided for MOT17-02-DPM

5. Pre-processing and post-processing is not part of this repository. We do provide pre-processed input detections for MOT17, though. 


## Evaluation

In order to run the tracker, call

`TrackingViaLiftedPath(sequenceName,datasetName,isTrain)`.

* `sequenceName`: name of the sequence.
* `datasetName`: name of the dataset of the sequence.
* `isTrain`: Boolean stating if the sequence is part of the training or test set.

A call can thus be:

`TrackingViaLiftedPath('MOT17-02-DPM','MOT17',true)`

## Graph Weights

Graph weight computations are not part of this repository. You can use your own definitions which must create a probability value for any edge within a defined time threshold (according to `LIFTED_SECONDS` and `BASE_SECONDS`). The function `computePairwiseCosts` should then load your pre-computed costs and call the solver on the resulting graph. The current format of the file containing pre-defined features can be seen in the `features/data/tmp/pairwise/fusion/` folder and is documented in the `computePairwiseCosts.m` function. 

## Input detections

Input detections are expected to be a CSV file with the column names "Frame,X_UL,Y_UL,W,H,PDET":

* `Frame`: the frame number of the detetection
* `X_UL`: x-coordinate of the detection box at the upper left corner
* `Y_UL`: y-coordinate of the detection box at the upper left corner
* `W`: width of the detection box
* `H`: height of the detection box
* `PDET`: probability of the detection box

