function trajectories = filter_trajectories(trajectories,paramsTracking,sequenceProperties)

rmIdx = [];
for iTraj = 1:numel(trajectories)
    
    for jTraj = iTraj+1:numel(trajectories)
       
        common_frames = intersect(trajectories(iTraj).descriptor_.fn,trajectories(jTraj).descriptor_.fn);

        if ~isempty(common_frames)
            iou_ratio = [];
            for fn = 1:numel(common_frames)
                [x_i,y_i,w_i,h_i] = trajectories(iTraj).selectDetectionsFromFrames(common_frames(fn)).getCoordinatesUpperLeft;
                [x_j,y_j,w_j,h_j] = trajectories(jTraj).selectDetectionsFromFrames(common_frames(fn)).getCoordinatesUpperLeft;
                iou_ratio(fn) = bboxOverlapRatio([x_i,y_i,w_i,h_i],[x_j,y_j,w_j,h_j]);
            end
            if mean(iou_ratio)>0.8
                % remove shorter trajectory if two trajectories have consistent high overlap, and if the shorter trajectory is completely covered 
                % by the larger trajectory
                if numel(trajectories(iTraj).descriptor_.fn) <= numel(trajectories(jTraj).descriptor_.fn) && numel(trajectories(iTraj).descriptor_.fn) == numel(common_frames) 
                    rmIdx = [rmIdx, iTraj];
                elseif numel(trajectories(jTraj).descriptor_.fn) <= numel(trajectories(iTraj).descriptor_.fn) && numel(trajectories(jTraj).descriptor_.fn) == numel(common_frames) 
                    rmIdx = [rmIdx, jTraj];
                end
            end
        end
    end
    
end
rmIdx = unique(rmIdx);
trajectories(rmIdx) = [];


end
