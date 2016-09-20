function [ad] = ucb(b, pHat, n_impressions, impressionsReceived, gamma)

score = b .* pHat + gamma * b .* sqrt(log(n_impressions) ./ impressionsReceived(1, :));
[~, ad] = max(score);

end

