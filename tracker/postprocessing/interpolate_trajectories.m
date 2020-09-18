function trajectories = interpolate_trajectories(trajectories,paramsTracking,sequenceProperties)


    for iTraj = 1:numel(trajectories)
        
        traj_fn = trajectories(iTraj).descriptor_.fn;
        start_fn = min(traj_fn);
        end_fn = max(traj_fn);
        interp_fn = start_fn:end_fn;
        
        ximg = trajectories(iTraj).descriptor_.ximg;
        yimg = trajectories(iTraj).descriptor_.yimg;
        width = trajectories(iTraj).descriptor_.width;
        height = trajectories(iTraj).descriptor_.height;
        pdet = trajectories(iTraj).descriptor_.pdet;
        identifier = trajectories(iTraj).descriptor_.identifier;
        ximg_interp = NaN(numel(interp_fn),1);
        yimg_interp = NaN(numel(interp_fn),1);
        width_interp = NaN(numel(interp_fn),1);
        height_interp = NaN(numel(interp_fn),1);
        pdet_interp = NaN(numel(interp_fn),1);
        
        
        ximg_interp(ismember(interp_fn,traj_fn)) = ximg;
        yimg_interp(ismember(interp_fn,traj_fn)) = yimg;
        width_interp(ismember(interp_fn,traj_fn)) = width;
        height_interp(ismember(interp_fn,traj_fn)) = height;
        pdet_interp(ismember(interp_fn,traj_fn)) = pdet;
        
        
        box_data = [ximg_interp,yimg_interp,width_interp,height_interp,pdet_interp];
        
        
        hasAssignment = ~isnan(ximg_interp);
        time_miss_assignments = interp_fn(~hasAssignment)/sequenceProperties.frameRate;
        time_assignments = interp_fn(hasAssignment)./sequenceProperties.frameRate;
        time_assignments_all = interp_fn./sequenceProperties.frameRate;
        
        startPos = 1;
        reset_idx = [];
        while ~isempty(find(hasAssignment(startPos:end) == 0,1))
            start_interp = (startPos-1)+find(hasAssignment(startPos:end) == 0,1)-1;
            end_interp = start_interp+find(hasAssignment(start_interp+1:end) == 1,1);
            startPos = end_interp;
            for component = 1:4
                box_data(start_interp+1:end_interp-1,component) = interp1(time_assignments_all([start_interp,end_interp]),box_data([start_interp,end_interp],component),time_assignments_all([start_interp+1:end_interp-1]),'linear');
            end
        end
        
        % replacing trajectory with interpolated trajectory
        box_data(hasAssignment,1) = trajectories(iTraj).descriptor_.ximg;
        box_data(hasAssignment,2) = trajectories(iTraj).descriptor_.yimg;
        box_data(hasAssignment,3) = trajectories(iTraj).descriptor_.width;
        box_data(hasAssignment,4) = trajectories(iTraj).descriptor_.height;
        box_data(isnan(box_data(:,5)),5) = -1;
        trajectories(iTraj).descriptor_.label = trajectories(iTraj).descriptor_.label(1)*ones(numel(interp_fn),1);
        trajectories(iTraj).descriptor_.fn = interp_fn';
        trajectories(iTraj).descriptor_.ximg = box_data(:,1);
        trajectories(iTraj).descriptor_.yimg = box_data(:,2);
        trajectories(iTraj).descriptor_.iou = 1*ones(numel(interp_fn),1);
        trajectories(iTraj).descriptor_.width = box_data(:,3);
        trajectories(iTraj).descriptor_.height = box_data(:,4);
        trajectories(iTraj).descriptor_.pdet= box_data(:,5);
        trajectories(iTraj).descriptor_.identifier= -1*ones(numel(interp_fn),1);
        trajectories(iTraj).descriptor_.identifier(hasAssignment) = identifier;
        trajectories(iTraj).descriptor_.interpolated(hasAssignment) = 1;
        
    end

end
