function TrackingViaLiftedPath(sequenceName,datasetName,isTrain)
%setup tracker
init_fol()

%folder name where configuration files are stored
dirConfigFiles = 'runtime_configuration';
[paramsTracking,configurations] = initTracker(dirConfigFiles);

% Register gurobi to Matlab
setup_gurobi(configurations)

% Load sequence properties and detections
mot_dataset = CMOT_Dataset.load_dataset(datasetName,isTrain,configurations);
mot_sequence = mot_dataset.load_sequence(sequenceName);


max_dist_time_base = paramsTracking.SOLVER.BASE_SECONDS;
paramsTracking.SOLVER.MAX_TIMEGAP_BASE = floor(max_dist_time_base*mot_sequence.get_sequence_properties().frameRate);

max_dist_time_lifted = paramsTracking.SOLVER.LIFTED_SECONDS;
paramsTracking.SOLVER.MAX_TIMEGAP_LIFTED = floor(max_dist_time_lifted*mot_sequence.get_sequence_properties().frameRate);

disp(['MAXIMAL TimeGap (Base) in Frames = ', num2str(paramsTracking.SOLVER.MAX_TIMEGAP_BASE)])
disp(['MAXIMAL TimeGap (Lifted) in Frames = ', num2str(paramsTracking.SOLVER.MAX_TIMEGAP_LIFTED)])

%%
% We store results dependent on paramters
% This allows to optimize paramter efficiently
% Add parameters which shall be analyzed

% Add here parameters which shall be studied that influence the solver
solver_study= ['ts',num2str(paramsTracking.SOLVER.TRACKLET_SIZE),'ab',num2str(paramsTracking.SOLVER.ALL_BASE_TRACKLET),'lthr',num2str(paramsTracking.SOLVER.POSITIVE_THRESHOLD_LIFTED)];

% Add here parameters which shall be studied that influence the graph
param_study = ['MAX_S-',num2str(paramsTracking.SOLVER.BASE_SECONDS),'IO',num2str(paramsTracking.SOLVER.INPUT_COST)];

solver_input_fol = fullfile(configurations.Tracker.outputDir,param_study,sequenceName);
paramsTracking.SOLVER.INPUT_GRAPH = fullfile(solver_input_fol,'problemDesc');
paramsTracking.SOLVER.INPUT_FRAMES = fullfile(solver_input_fol,'problemDesc_frames');
paramsTracking.SOLVER.OUTPUT_PATH = fullfile(solver_input_fol,solver_study,'/');

create_folder(paramsTracking.SOLVER.OUTPUT_PATH);


% the solver is run on a copy of configuration file
configFile = tempname;
write_params_to_config(paramsTracking,configFile);
paramsTracking.SOLVER.filename = configFile;


saveDir = paramsTracking.SOLVER.OUTPUT_PATH ;

disp(['Storing results into ',saveDir])


saveImg = configurations.Tracker.saveImg;

if saveImg
    visObj = CVisualizer(mot_sequence,saveImg,saveDir);
else
    visObj = [];
end

%%
saveDir_solver = fullfile(saveDir,'solver');
if exist(saveDir_solver,'dir') ~= 7
    mkdir(saveDir_solver)
end

configurations.saveDir_solver = saveDir_solver;
%%

weightedNodes = computeUnaryCosts(mot_sequence,paramsTracking);
if configurations.Solver.CREATE_PROBLEM_FILE
    
    weightedEdgesBase = computePairwiseCosts(mot_sequence,paramsTracking);
    % possible to implement function specifically for lifted edges
    %weightedEdgesLifted = computeLiftedCosts(detections,paramsTracking,sequenceProperties);
    weightedEdgesLifted = [];
    
    
    problemGraph.nodes = weightedNodes;
    problemGraph.base_edges = weightedEdgesBase;
    problemGraph.lifted_edges = weightedEdgesLifted;
end
associations = solveLiftedDP(problemGraph,mot_sequence,paramsTracking);
disp('solver done')

if configurations.Solver.EVAL_SOLVER
    clear LDPtrajectories;
    detections = mot_sequence.detections_;
    for iTraj = 1:numel(associations)
        nodeList = weightedNodes.convert_graphNodeID_2_boxID(associations{iTraj});
        LDPtrajectories(iTraj) = detections.selectDetectionsFromBoxID(nodeList);
        LDPtrajectories(iTraj).descriptor_.label = iTraj*ones(numel([LDPtrajectories(iTraj).descriptor_.label]),1);
    end
    disp('Saving results into mat file')
    save(fullfile(saveDir,'trajectories.mat'),'LDPtrajectories','associations','visObj','mot_sequence','problemGraph','paramsTracking','saveDir','-v7.3');
    disp('Saving done')
    
    if saveImg
        disp('Creating association video')
        visObj.drawTrajectories(LDPtrajectories,detections);
        [status, result] = system(['ffmpeg -y -i ',saveDir,'%06d.png ',saveDir,'/',sequenceName,'_association.mp4']);
        [status, result] = system(['rm  ',saveDir,'*.png ']);
        disp('association video done')
    end
    
    saveDir_interp = fullfile(saveDir,'interp/');
    if exist(saveDir_interp,'dir') ~= 7
        mkdir(saveDir_interp)
    end
    
    configurations.saveDir = saveDir_interp;
    
    disp('Linear interpolation post processing')
    trajectoriesAfterPostProcessing = postProcessing(LDPtrajectories,paramsTracking,mot_sequence);
    disp('Linear interpolation post processing done ')
    result_file = fullfile(saveDir_interp,mot_sequence.get_sequence_properties().name);
    CNodes.saveAsCSV(trajectoriesAfterPostProcessing,result_file)
 
    if saveImg
        visObj = CVisualizer(mot_sequence,saveImg,saveDir_interp);
        visObj.drawTrajectories(trajectoriesAfterPostProcessing,detections);
        [status, result] = system(['ffmpeg -y -i ',saveDir_interp,'%06d.png ',saveDir_interp,sequenceName,'.mp4']);
        [status, result] = system(['rm  ',saveDir_interp,'*.png ']);
    end
    
end
end
