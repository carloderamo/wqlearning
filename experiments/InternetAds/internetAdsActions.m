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
n_actions = 10:10:100;
error = zeros(n_experiments, 3, length(n_actions));

bias = zeros(length(n_actions), 3);
variance = zeros(length(n_actions), 3);
mse = zeros(length(n_actions), 3);

for i = 1:length(n_actions)
    n_trials = 10000;
    
    idxs = repmat(1:n_actions(i), n_actions(i), 1);
    idxs = (1 - eye(n_actions(i))) .* idxs';
    idxs(idxs == 0) = [];
    idxs = reshape(idxs, n_actions(i) - 1, n_actions(i));
    idxs = idxs';

    for experiment = 1:n_experiments

        fprintf('Experiment: %d\n', experiment);

	p = 0.02 + (0.05 - 0.02) * rand(1, n_actions(i));

        [clicks, means, sigma] = crtLearning(n_actions(i), n_trials, p);

        % Maximum Estimator
        error(experiment, 1, i) = max(sum(clicks) / n_trials) - max(p);

        % Double Estimator
        clicks1 = clicks(1:floor(size(clicks, 1) / 2), :);
        clicks2 = clicks(floor(size(clicks, 1) / 2) + 1:end, :);
        pHat1 = sum(clicks1) / size(clicks1, 1);
        pHat2 = sum(clicks2) / size(clicks2, 1);
        doubleAdMax = find(pHat1 == max(pHat1));
        mu1 = mean(pHat2(doubleAdMax));

        doubleAdMax = find(pHat2 == max(pHat2));
        mu2 = mean(pHat1(doubleAdMax));
        error(experiment, 2, i) = mean([mu1 mu2]) - max(p);

        % W Estimator
        %lower_limit = means - 8 * sigma;
        %upper_limit = means + 8 * sigma;
        %n_trapz = 1e2;
        %x = zeros(n_trapz, n_actions(i));
        %y = zeros(size(x));
        %for j = 1:n_actions(i)
        %    x(:, j) = linspace(lower_limit(j), upper_limit(j), n_trapz);
        %    y(:, j) = normpdf(x(:, j), means(j), sigma(j)) .* ...
        %            prod(normcdf(repmat(x(:, j), 1, n_actions(i) - 1), ...
        %                         means(repmat(idxs(j, :), n_trapz, 1)), ...
        %                         sigma(repmat(idxs(j, :), n_trapz, 1))), 2);
        %end
        %integrals = trapz(y, 1) .* ((upper_limit - lower_limit) / (n_trapz - 1));
        %error(experiment, 3, i) = integrals * means' - max(p);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% SAMPLING PROBS %%%%%%
	n_samples = 1000;
	samples = repmat(means, n_samples, 1) + repmat(sigma, n_samples, 1) .* randn(n_samples, n_actions(i));
	[~, max_idxs] = max(samples');
	max_count = zeros(size(samples));
	max_count(sub2ind(size(samples), 1:length(max_idxs'), max_idxs)) = 1;

	probs = sum(max_count, 1) / n_samples;
	error(experiment, 3, i) = probs * means' - max(p);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    bias(i, :) = mean(error(:, :, i));
    variance(i, :) = var(error(:, :, i));
    mse(i, :) = bias(i, :).^2 + variance(i, :);
end

n_actions_text = 'actions';
save(strcat('internetAds-', n_actions_text,'.mat'));
