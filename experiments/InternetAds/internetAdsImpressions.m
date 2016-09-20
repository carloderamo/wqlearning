%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script aims to replicate the internet ads   %
% experiments presented in H. V. Hasselt's article %
% "Estimating the Maximum Expected Value: An       %
% Analysis of (Nested) Cross Validation and the    %
% Maximum Sample Average".                         %
%                                                  %
% Written by: Carlo D'Eramo                        %
%                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

n_experiments = 2000;
n_actions = 30;
n_visitors = 30000:30000:300000;
error = zeros(n_experiments, 3, length(n_visitors));

bias = zeros(length(n_visitors), 3);
variance = zeros(length(n_visitors), 3);
mse = zeros(length(n_visitors), 3);

idxs = repmat(1:n_actions, n_actions, 1);
idxs = (1 - eye(n_actions)) .* idxs';
idxs(idxs == 0) = [];
idxs = reshape(idxs, n_actions - 1, n_actions);
idxs = idxs';

for i = 1:length(n_visitors)
    n_trials = floor(n_visitors(i) / n_actions);
    
    for experiment = 1:n_experiments

        fprintf('Experiment: %d\n', experiment);

	p = 0.02 + (0.05 - 0.02) * rand(1, n_actions);

        [clicks, means, sigma] = crtLearning(n_actions, n_trials, p);

        % Maximum Estimator
        error(experiment, 1, i) = max(sum(clicks) / n_trials) - max(p);

        % Double Estimator
        clicks1 = clicks(1:floor(n_trials / 2), :);
        clicks2 = clicks(floor(n_trials / 2) + 1:end, :);
        pHat1 = sum(clicks1) / floor(n_trials / 2);
        pHat2 = sum(clicks2) / floor(n_trials / 2);
        doubleAdMax = find(pHat1 == max(pHat1));
        mu1 = mean(pHat2(doubleAdMax));

        doubleAdMax = find(pHat2 == max(pHat2));
        mu2 = mean(pHat1(doubleAdMax));
        error(experiment, 2, i) = mean([mu1 mu2]) - max(p);

        % W Estimator
        lower_limit = means - 8 * sigma;
        upper_limit = means + 8 * sigma;
        n_trapz = 1e2;
        x = zeros(n_trapz, n_actions);
        y = zeros(size(x));
        for j = 1:n_actions
            x(:, j) = linspace(lower_limit(j), upper_limit(j), n_trapz);
            y(:, j) = normpdf(x(:, j), means(j), sigma(j)) .* ...
                    prod(normcdf(repmat(x(:, j), 1, n_actions - 1), ...
                                 means(repmat(idxs(j, :), n_trapz, 1)), ...
                                 sigma(repmat(idxs(j, :), n_trapz, 1))), 2);
        end
        integrals = trapz(y, 1) .* ((upper_limit - lower_limit) / (n_trapz - 1));
        error(experiment, 3, i) = integrals * means' - max(p);
    end
    
    bias(i, :) = mean(error(:, :, i));
    variance(i, :) = var(error(:, :, i));
    mse(i, :) = bias(i, :).^2 + variance(i, :);
end

n_actions_text = 'impressions';
save(strcat('internetAds-', n_actions_text,'.mat'));

