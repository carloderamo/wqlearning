function [pHat, clicksHistory, means, sigma] = debiasedCRT(p, b, n_impressions, gamma, adaptive)

sumPTimesB = zeros(1, length(p));
sumSquarePTimesB = zeros(1, length(p));

pHat = zeros(1, length(p));
clicks = zeros(2, length(p));
clicksHistory = zeros(floor(n_impressions / 2), length(p), 2);
impressionsReceived = zeros(1, length(p));
for i = 1:floor(n_impressions / 2)
    if adaptive
        ad = ucb(b, pHat, n_impressions, impressionsReceived, gamma);
    else
        ad = uniform(b, pHat, gamma);
    end
    
    click_event = 0; % used for W estimator
    impressionsReceived(1, ad) = impressionsReceived(1, ad) + 1;
    clicksHistory(i, ad, 2) = 2; % used for count impressions in the double estimator case in setting 4
    for j = 1:2
        if rand < p(ad)
            clicks(j, ad) = clicks(j, ad) + 1;
            clicksHistory(i, ad, 1) = clicksHistory(i, ad, 1) + 1;
            click_event = click_event + 1;
        end
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

pHat = clicks(2, :) ./ impressionsReceived;
pHat(isnan(pHat)) = 0;

end