classdef CWeightedEdges < handle
    
    properties
        costs_
        sourceIdentifier_
        sinkIdentifier_
        sourceFrame_
        sinkFrame_
    end
    
    
    methods (Static = true)
        
        function weightedEdges = weightedEdgesFromMatrix(weightedMat,identifiersRows,identifiersCols)
            costs = [];
            if issparse(weightedMat)
                [rMat,cMat] = find(weightedMat ~= 0);
            else
                [rMat,cMat] = find(~isnan(weightedMat));
            end
            for iEntry = 1:numel(rMat)
                costs = [costs, weightedMat(rMat(iEntry),cMat(iEntry))];
                nodeIdentSource = [nodeIdentSource, identifiersRows(rMat(iEntry))];
                nodeIdentSink = [nodeIdentSink, identifiersCols(cMat(iEntry))];
            end
            weightedEdges = CWeightedEdges(nodeIdentSource,nodeIdentSink,costs); 
        end     
    end
    
    methods
        
        function obj = CWeightedEdges(nodeIdentifiersSource,nodeIdentifiersSink,costs,fnSource,fnSink)
            assert(numel(nodeIdentifiersSink) == numel(nodeIdentifiersSource) && numel(nodeIdentifiersSink) == numel(costs),' Number of source, sink and costs for edge created must match');
            obj.costs_ = costs;
            obj.sourceIdentifier_ = nodeIdentifiersSource;
            obj.sinkIdentifier_ = nodeIdentifiersSink;
            obj.sourceFrame_ = fnSource;
            obj.sinkFrame_ = fnSink;
        end
        
        
        function [nodeIDSource,nodeIDSink,edgeCosts,fnSource,fnSink] = getEdgeDescriptors(obj,nodeIdx)
            if nargin == 1
                nodeIdx = 1:numel(obj.costs_);
            end
            nodeIDSource = obj.sourceIdentifier_(nodeIdx);
            nodeIDSink = obj.sinkIdentifier_(nodeIdx);
            edgeCosts = obj.costs_(nodeIdx);
            fnSource = obj.sourceFrame_(nodeIdx);
            fnSink = obj.sinkFrame_(nodeIdx);
        end 
        
        
        
    end
end
