function [pHat, means, sigma] = biasedCRT(p, b, n_impressions, gamma, adaptive)

sumPTimesB = zeros(1, length(p));
sumSquarePTimesB = zeros(1, length(p));

pHat = zeros(1, length(p));
clicks = zeros(1, length(p));
impressionsReceived = zeros(1, length(p));

for i = 1:n_impressions
    if adaptive
        ad = ucb(b, pHat, n_impressions, impressionsReceived, gamma);
    else
        ad = uniform(b, pHat, gamma);
    end    
    
    click_event = 0; % used for W estimator
    impressionsReceived(1, ad) = impressionsReceived(1, ad) + 1;
    if rand < p(ad)
        clicks(1, ad) = clicks(1, ad) + 1;
        click_event = 1;
    end
    pHat(1, ad) = clicks(1, ad) / impressionsReceived(1, ad);

    sumPTimesB(1, ad) = sumPTimesB(1, ad) + b(1, ad) * click_event;
    sumSquarePTimesB(1, ad) = sumSquarePTimesB(1, ad) + (b(1, ad) * click_event)^2;
end

means = sumPTimesB ./ impressionsReceived;
means(isnan(means)) = 0;
variance = (sumSquarePTimesB ./ impressionsReceived - means.^2 ) ./ impressionsReceived;
variance(variance < 0) = 0;
sigma = sqrt(variance);
sigma(isnan(sigma)) = 1e10;
sigma(sigma == 0) = 1e-4;

pHat = clicks ./ impressionsReceived;
pHat(isnan(pHat)) = 0;

end