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

n_experiments = 5000;
n_actions = 30;
n_visitors = 3e5;
n_trials = n_visitors / n_actions;

idxs = repmat(1:n_actions, n_actions, 1);
idxs = (1 - eye(n_actions)) .* idxs';
idxs(idxs == 0) = [];
idxs = reshape(idxs, n_actions - 1, n_actions);
idxs = idxs';

max_p = 0.02:0.01:0.1;

error = zeros(n_experiments, 4, length(max_p));
bias = zeros(length(max_p), 4);
variance = zeros(length(max_p), 4);
mse = zeros(length(max_p), 4);

for i = 1:length(max_p)
    errorME = zeros(n_experiments, 1);
    errorDE = zeros(n_experiments, 1);
    errorMME = zeros(n_experiments, 1);
    errorWE = zeros(n_experiments, 1);
    parfor experiment = 1:n_experiments

        fprintf('N. setting: %d, Experiment: %d\n', i, experiment);

	    p = 0.02 + (max_p(i) - 0.02) * rand(1, n_actions);

        [clicks, means, sigma] = crtLearning(n_actions, n_trials, p);

        % Maximum Estimator
        errorME(experiment, 1) = max(sum(clicks) / n_trials) - max(p);

        % Double Estimator
        clicks1 = clicks(1:n_trials / 2, :);
        clicks2 = clicks(n_trials / 2 + 1:end, :);
        pHat1 = sum(clicks1) / (n_trials / 2);
        pHat2 = sum(clicks2) / (n_trials / 2);
        doubleAdMax = find(pHat1 == max(pHat1));
        mu1 = mean(pHat2(doubleAdMax));

        doubleAdMax = find(pHat2 == max(pHat2));
        mu2 = mean(pHat1(doubleAdMax));
        errorDE(experiment, 1) = mean([mu1 mu2]) - max(p);
        
         % Maxmin Estimator
        clicks1 = clicks(1:n_trials / 2, :);
        clicks2 = clicks(n_trials / 2 + 1:end, :);
        mu1 = sum(clicks1) / (n_trials / 2);
        mu2 = sum(clicks2) / (n_trials / 2);
        mumin = min(mu1, mu2);
        errorMME(experiment, 1) = max(mumin) - max(p);

        % W Estimator
        lower_limit = means - 8 * sigma;
        upper_limit = means + 8 * sigma;
        n_trapz = 3e2;
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
        errorWE(experiment, 1) = integrals * means' - max(p);
    end
    error = cat(2, errorME, errorDE, errorMME, errorWE);
    bias(i, :) = mean(error);
    variance(i, :) = var(error);
    mse(i, :) = bias(i, :).^2 + variance(i, :);
end

n_actions_text = 'probs';
save(strcat('internetAds-', n_actions_text,'.mat'));
