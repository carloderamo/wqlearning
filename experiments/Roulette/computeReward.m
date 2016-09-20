function [reward] = computeReward(actionN)
% This function sample the winning event with the probability
% related to the action.

reward = 0;

if actionN >= 1 && actionN <= 38
    nSquares = 1;
    reward = sampleOutcome(nSquares);
elseif actionN >= 39 && actionN <= 96
    nSquares = 2;
    reward = sampleOutcome(nSquares);
elseif actionN >= 97 && actionN <= 111
    nSquares = 3;
    reward = sampleOutcome(nSquares);
elseif actionN >= 112 && actionN <= 133
    nSquares = 4;
    reward = sampleOutcome(nSquares);
elseif actionN == 134
    nSquares = 5;
    reward = sampleOutcome(nSquares);
elseif actionN >= 135 && actionN <= 145
    nSquares = 6;
    reward = sampleOutcome(nSquares);
elseif actionN >= 146 && actionN <= 151
    nSquares = 12;
    reward = sampleOutcome(nSquares);
elseif actionN >= 152 && actionN <= 157
    nSquares = 18;
    reward = sampleOutcome(nSquares);

end

end

