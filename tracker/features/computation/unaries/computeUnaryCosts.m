function nodes = computeUnaryCosts(mot_sequence,paramsTracking)

detections = mot_sequence.detections_;
costsConf = paramsTracking.SOLVER.CONF_COST*[detections.descriptor_.pdet];

costsIn = paramsTracking.SOLVER.INPUT_COST*ones(size(costsConf,1),size(costsConf,2));
costsOut = paramsTracking.SOLVER.OUTPUT_COST*ones(size(costsConf,1),size(costsConf,2));

nodes = CWeightedNodes(detections,costsConf,costsIn,costsOut);
end