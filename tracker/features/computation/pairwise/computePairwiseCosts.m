function weightedEdges = computePairwiseCosts(mot_sequence,paramsTracking)
% fusion costs format:
% CSV table with columns 
% prob_f : probabilty of an edge (1 = connect nodes, 0 = disconnect nodes)
% pw_info_0, pw_info_1 : node indices of an edge, indices start with 0, 
% each index+1 denotes the row number of the detection file corresponding
% to the node
% pw_info_2 - pw_info_4 : they can be ignored
%
% Soft constraints may be used, which punish connection of nodes
% soft constrain format:
% CSV table with columns:
% reset_edges_0, reset_edges_1 : node indices of an edge, indices start
% with 0. Each index+1 denotes the row number of the detection file corresponding
% to the node

softConstraints = true;


detections = mot_sequence.detections_;
sequenceProperties = mot_sequence.get_sequence_properties();
configurations = mot_sequence.get_configurations();
[source,sink] = linksWithinWindow(detections,paramsTracking.SOLVER.MAX_TIMEGAP_BASE);

edge_list = [source,sink];
edge_list = sortrows(edge_list);
source = edge_list(:,1);
sink = edge_list(:,2);

identities = detections.getIdentifiers();
identifierSource = identities(source);
identifierSink= identities(sink);

fn = detections.frameRange;

fusion_file = fullfile(configurations.Data.featuresDir,'pairwise/fusion/',['fusion_costs_IN=',sequenceProperties.name,'.csv']);
pw_file = dir(fusion_file);
assert(numel(pw_file) == 1)
pairwise_data = readtable(fusion_file);



if softConstraints
    soft_file = fullfile(configurations.Data.featuresDir,'pairwise/fusion/',['reset_fusion_costs_IN=',sequenceProperties.name,'.csv']);
    reset_file = dir(soft_file);
    assert(numel(reset_file) == 1)
    reset_edges = table2array(readtable(soft_file));
    reset_edges = reset_edges+1;
end
edge_data = [pairwise_data.pw_info_0,pairwise_data.pw_info_1];

edge_data = edge_data(:,1:2)+1;
edge_costs = prob2logit(pairwise_data.prob_f);

if softConstraints
    pairw_costs = [edge_costs; repmat(15,size(reset_edges,1),1)];
    edge_info_fusion = [edge_data; reset_edges];
else
    edge_info_fusion = edge_data;
end

[~,index_list_to_intersect,index_features_to_intersect] = intersect([identifierSource, identifierSink],edge_info_fusion,'rows');
assert(size(index_list_to_intersect,1) == size(edge_list,1),'features missing for some edges')
pairw_costs = double(pairw_costs(index_features_to_intersect));

weightedEdges = CWeightedEdges(identifierSource,identifierSink,pairw_costs,fn(source),fn(sink));
end
