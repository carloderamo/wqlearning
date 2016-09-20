function [ad] = uniform(b, pHat, eps)

score = b .* pHat;

if rand < 1 - eps
    ad = datasample(find(score == max(score)), 1);
else
    ad = floor(rand * length(b)) + 1;

end

