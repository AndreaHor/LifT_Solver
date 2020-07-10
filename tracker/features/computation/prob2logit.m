function costs = prob2logit(probs)
%PROB2LOGIT Summary of this function goes here
%   Detailed explanation goes here
costs = log((1-probs)./probs);
costs(costs<-10) = -10;
costs(costs>100) = 100;
end

