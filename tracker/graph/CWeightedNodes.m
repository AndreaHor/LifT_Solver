classdef CWeightedNodes < handle
    % Each node consists of
    % boxID_ : unique identifier of each detection, starting from 1,
    % which is the row number of the detection
    % costsConf_ : confidence costs of detection
    % costsIn_   : in costs (start of a trajectory at the detection)
    % costsOut_ : out costs (end of a trajectory at the detection)
    %
    % Each boxID_ has a corresponding graphID value (one-to-one)
    % graphID is ensured to be consecutive numbers starting from 0.
    % It is used in the graph description file of the lifted disjoint paths
    % solver.
    % Note: If a subselection on the detections is performed boxID_ might not
    % consective numbers
    %
    % The mapping between boxID_ and graphID is given via
    % boxID2graphNodeID_ : maps uid_ to graphID
    % graphNodeID2boxID_ : maps graphID to uid_
    
    properties (Access = private)
        boxID_
        costsConf_
        costsIn_
        costsOut_
        graphNodeID2boxID_
        boxID2graphNodeID_
        frames_
    end
    
    methods
        function obj = CWeightedNodes(detections,costsConf,costsIn,costsOut)
            uids = detections.getIdentifiers();
            [fn,~,~] = detections.frameRange;
            obj.frames_ = fn;
            obj.boxID_ = [[uids(:)]'];
            obj.costsConf_ = [[costsConf(:)]'];
            if ~isempty(costsIn) && ~isempty(costsOut)
                obj.costsIn_ = [costsIn];
                obj.costsOut_ = [costsOut];
            else
                obj.costsIn_ = [NaN(1,numel(uids))];
                obj.costsOut_ = [NaN(1,numel(uids))];
            end
            
            %initialize mappings
            obj.get_node_descriptors();
        end
        
        function nNodes = getNumberOfNodes(obj)
            nNodes = numel(obj.boxID_);
        end
        
        function nodeIds = convert_boxID_2_graphNodeID(obj,boxIDVec)
            nodeIds = obj.boxID2graphNodeID_(boxIDVec)-1;
            assert(~any(isnan(nodeIds)),'Some boxes are not part of this node object.')
        end
        
        function boxIds = convert_graphNodeID_2_boxID(obj,nodeVec)
            boxIds = obj.graphNodeID2boxID_(nodeVec+1);
            assert(~any(isnan(boxIds)),'Some nodes are missing corresponding detection boxes.')
        end
        
        function [boxID,nodeCost,nodeStartCosts,nodeEndCosts,frames] = get_node_descriptors(obj,box_selection)
            % box_selection: set of detection indices to get node
            % descriptor for
            if nargin == 1
                box_selection = 1:numel(obj.boxID_);
            end
            boxID    = obj.boxID_(box_selection);
            nodeCost        = obj.costsConf_(box_selection);
            nodeStartCosts  = obj.costsIn_(box_selection);
            nodeEndCosts  = obj.costsOut_(box_selection);
            
            [~,~,nodeID] = unique(boxID);
            
            assert(issorted(boxID) && numel(boxID) == numel(unique(boxID)),'Node object has some detections multiple times')
            obj.graphNodeID2boxID_ = NaN(max(nodeID),1);
            obj.graphNodeID2boxID_(nodeID(:)) = boxID(:);
            
            obj.boxID2graphNodeID_= NaN(max(boxID),1);
            obj.boxID2graphNodeID_(boxID(:)) = nodeID(:);
            
            frames = obj.frames_;
        end
    end
end