function associationList = solveLiftedDP(graph,mot_sequence,paramsTracking)


configurations = mot_sequence.get_configurations();

solverIODir = [paramsTracking.SOLVER.OUTPUT_PATH];

if configurations.Solver.CREATE_PROBLEM_FILE
    weightedNodes = graph.nodes;
    weightedEdgesNormal = graph.base_edges;
    weightedEdgesLifted = graph.lifted_edges;
    
    
    %CLEAN FILES
    delete([solverIODir,'/',paramsTracking.SOLVER.OUTPUT_PREFIX,'*'])
    delete(paramsTracking.SOLVER.INPUT_GRAPH)
    delete(paramsTracking.SOLVER.INPUT_FRAMES)
    
    %Creating graph file
    inputFileFullPath = fullfile(paramsTracking.SOLVER.INPUT_GRAPH);
    fid = fopen(inputFileFullPath,'w');
    if fid == -1
        error(['Not able to create file: ',inputFileFullPath])
    end
    % adding unaries to problem description
    [boxID,nodeCost,nodeStartCosts,nodeEndCosts,nodeFrames] = weightedNodes.get_node_descriptors();
    nodeID = weightedNodes.convert_boxID_2_graphNodeID(boxID);
    
    
    disp('writing nodes to input file')
    fprintf(fid,'%d \n',numel(nodeID));
    for iNode = 1:numel(nodeID)
        % setting in / out costs via solver params
        if abs(nodeCost(iNode))>0
            fprintf(fid,'%d,%f \n',nodeID(iNode),nodeCost(iNode));
        end
    end
    fprintf(fid,'\n');
    
    % adding pairwise costs to problem description TODO VERIFY BOX ID and NODE ID is used correctly
    
    disp('writing edges to input file')
    [boxIDSource,boxIDSink,edgeCosts,fnSource,fnSink] = weightedEdgesNormal.getEdgeDescriptors();
    nodeIDStart = weightedNodes.convert_boxID_2_graphNodeID(boxIDSource);
    nodeIDEnd   = weightedNodes.convert_boxID_2_graphNodeID(boxIDSink);
    for iEdge = 1:numel(nodeIDStart)
        fprintf(fid,'%d,%d,%f \n',nodeIDStart(iEdge),nodeIDEnd(iEdge),edgeCosts(iEdge));
    end
    
    % adding lifted costs to problem description
    if ~isempty(weightedEdgesLifted)
        fprintf(fid,'\n');
        disp('writing lifted edges to input file')
        [nodeIDStart,nodeIDEnd,edgeCosts] = weightedEdgesLifted.getEdgeDescriptors();
        nodeIDStart = weightedNodes.convert_boxID_2_graphNodeID(nodeIDStart);
        nodeIDEnd   = weightedNodes.convert_boxID_2_graphNodeID(nodeIDEnd);
        for iEdge = 1:numel(nodeIDStart)
            fprintf(fid,'%d,%d,%f \n',nodeIDStart(iEdge),nodeIDEnd(iEdge),edgeCosts(iEdge));
        end
    end
    fclose(fid);
    disp('Input files created')
    
    timeFileFullPath = fullfile(paramsTracking.SOLVER.INPUT_FRAMES);
    
    
    disp('Creating frames file')
    fid = fopen(timeFileFullPath,'w');
    for iNode = 1:numel(nodeID)
        fprintf(fid,'%d,%d \n',nodeID(iNode),nodeFrames(iNode));
    end
    fclose(fid);
    
    solverDir = configurations.Solver.solverDir;
    solverName = 'run-disjont-paths';
    disp('calling lifted disjoint paths solver')
end

movefile(paramsTracking.SOLVER.filename,fullfile(solverIODir,'params_sequence.ini'))

if configurations.Solver.CALL_SOLVER
    
    solverCallStr = [fullfile(solverDir,solverName),' -c ',fullfile(solverIODir,'params_sequence.ini')];
    
    disp(['calling Solver via: ',solverCallStr])
    tic
    [status, result] = system(solverCallStr);
    toc
end

if configurations.Solver.EVAL_SOLVER
    searchKey = [paramsTracking.SOLVER.OUTPUT_PREFIX,'-all-paths-FINAL'];
    disp([solverIODir,searchKey,'.txt'])
    
    searchStr = fullfile(solverIODir,[searchKey,'.txt']);
    resultName =  dir(searchStr);
    if length(resultName) == 0
       result
       error('No output file found');
    end
    solutionMat = dlmread([fullfile(solverIODir,resultName.name)]);
    nextTraj = 0;
    
    for iRow = 1:size(solutionMat,1)
        nextTraj = nextTraj+1;
        endEntries = find(solutionMat(iRow,2:end) == 0,1);
        if ~isempty(endEntries)
            endcol = endEntries;
        else
            endcol = size(solutionMat,2);
        end
        associationList{nextTraj} = solutionMat(iRow,1:endcol);
    end
    
    disp('Zipping solver files')
    [~,~] = system(['cd ',solverIODir,'  && zip solver_data.zip',' *']);
    disp('Zipping solver files done')
    
    disp('Cleaning temp solver files')
    delete([solverIODir,'/',paramsTracking.SOLVER.OUTPUT_PREFIX,'*'])
    disp('Cleaning temp solver files done')

    movefile(fullfile(solverIODir,'solver_data.zip'),fullfile(solverIODir,'solver','solver_data.zip'));
else
    associationList = 0;
end
end

