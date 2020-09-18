function trajectories = postProcessing(trajectories,weightedEdgesBase,paramsTracking,sequenceProperties,configurations)

disp('interpolate trajectories')
trajectories  = interpolate_trajectories(trajectories,paramsTracking,sequenceProperties);
disp('interpolate trajectories done')

disp('filter trajectories')
trajectories = filter_trajectories(trajectories,paramsTracking,sequenceProperties);
disp('filter trajectories done')




end
