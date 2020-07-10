function [source,sink] = linksWithinWindow(detections,maxTimeGap)

frames = detections.frameRange;

uFrames = unique(frames);

linear_idx = 0;
for iFrame = 1:numel(uFrames)-1
    idx_current_nodes = find(frames == uFrames(iFrame)); %source nodes
    for jFrame = (iFrame+1):numel(uFrames)
        if uFrames(jFrame)-uFrames(iFrame) <= maxTimeGap && uFrames(jFrame)-uFrames(iFrame) > 0
            linear_idx = linear_idx+1;
            idx_future_nodes = find(frames == uFrames(jFrame)); % sink nodes        
            edge_list = combvec(idx_future_nodes',idx_current_nodes');
            sink_tmp{linear_idx} = edge_list(1,:);
            source_tmp{linear_idx} = edge_list(2,:);
        end
    end
    
end
sink = [sink_tmp{:}];
sink = sink';
source = [source_tmp{:}];
source = source';
end