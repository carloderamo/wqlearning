function [outcome] = sampleOutcome(nSquares)
% Sample the win event  with the probability
% computed w.r.t. the number of squares occupied
% by the bet.

bet = 1;
p = nSquares / 38;
payout = 36 / nSquares - 1;

if(rand('double') < p)
    outcome = bet * payout;
else
    outcome = -bet;
end

end

