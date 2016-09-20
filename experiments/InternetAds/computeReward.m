function [reward] = computeReward(p)
% This function sample the click event with the probability
% related to the action.

reward = 0;
if(rand < p)
    reward = 1;

end

