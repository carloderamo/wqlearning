%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script aims to replicate the sponsored search %
% experiments presented in Xu et al. article         %
% "Estimation Bias in Multi-Armed Bandit Algorithms  %
% for Search Advertising".                           %
%                                                    % 
% Written by: Carlo D'Eramo                          %
%                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

n_experiments = 2000;
gamma = 1;
settings = [1, 2, 3, 4];

%Iterate over the selected settings
for setting = 1:length(settings)
    whichSetting = settings(setting);
    
    
    if whichSetting ~= 3 && whichSetting ~= 4
        if whichSetting == 1
            p = [0.09 0.1]; b = [2 1]; n_impressions = [0.5 1 3 5 10 15] * 1e3;
        elseif whichSetting == 2
            p = [0.3 0.1]; b = [0.7 1]; n_impressions = [0.5 1 3 5 10 15] * 1e3;
        end
    	
	% Specific for Weighted estimator    
        idxs = repmat(1:length(p), length(p), 1);
        idxs = (1 - eye(length(p))) .* idxs';
        idxs(idxs == 0) = [];
        idxs = reshape(idxs, length(p) - 1, length(p));
        idxs = idxs';

        revenueNoDebias = zeros(n_experiments, length(n_impressions));
        revenueDebias = zeros(n_experiments, length(n_impressions));
        revenueWDebias = zeros(n_experiments, length(n_impressions));
        for experiment = 1:n_experiments
            fprintf('Experiment: %d\n', experiment);

            for i = 1:length(n_impressions)
                % Biased CTR Learning Phase
                [pHat, means, sigma] = biasedCRT(p, b, n_impressions(i), gamma, true);

                % SP Auction Phase
                pHatTimesB = b .* pHat;
                maxHatIdx = datasample(find(pHatTimesB == max(pHatTimesB)), 1);
                secondMaxHatIdx = datasample(find(pHatTimesB == max(pHatTimesB([1:maxHatIdx - 1, maxHatIdx + 1:end]))), 1);

                if pHat(maxHatIdx) == 0
                    revenueNoDebias(experiment, i) = 0;
                else
                    revenueNoDebias(experiment, i) = p(maxHatIdx) * b(secondMaxHatIdx) * pHat(secondMaxHatIdx) / pHat(maxHatIdx);
                end

                % Debiased CTR Learning Phase
                [pHatDouble, ~, ~, ~] = debiasedCRT(p, b, n_impressions(i), gamma, true);

                % SP Auction Phase
                pHatTimesB = b .* pHatDouble;
                maxHatIdx = datasample(find(pHatTimesB == max(pHatTimesB)), 1);
                secondMaxHatIdx = datasample(find(pHatTimesB == max(pHatTimesB([1:maxHatIdx - 1, maxHatIdx + 1:end]))), 1);

                if pHatDouble(maxHatIdx) == 0
                    revenueDebias(experiment, i) = 0;
                else
                    revenueDebias(experiment, i) = p(maxHatIdx) * b(secondMaxHatIdx) ...
                                                   * pHatDouble(secondMaxHatIdx) / pHatDouble(maxHatIdx);
                end
                
                % Weighted Debiased CTR Learning Phase
                % Same as no selection debiasing!

                % SP Auction Phase
                lower_limit = means - 8 * sigma;
                upper_limit = means + 8 * sigma;
                n_trapz = 1e2;
                x = zeros(n_trapz, length(p));
                y = zeros(size(x));
                
                for j = 1:length(p)
                    x(:, j) = linspace(lower_limit(j), upper_limit(j), n_trapz);
                    y(:, j) = normpdf(x(:, j), means(j), sigma(j)) .* ...
                            prod(normcdf(repmat(x(:, j), 1, length(p) - 1), ...
                                         means(repmat(idxs(j), n_trapz, 1))', ...
                                         sigma(repmat(idxs(j), n_trapz, 1))'), 2);
                end
                integralsMaximum = trapz(y, 1) .* ((upper_limit - lower_limit) / (n_trapz - 1));
                
                x = zeros(n_trapz, length(p));
                y = zeros(size(x));
                for j = 1:length(p)
                    x(:, j) = linspace(lower_limit(j), upper_limit(j), n_trapz);
                    cdfs = normcdf(repmat(x(:, j), 1, length(p) - 1), ...
                                   means(repmat(idxs(j), n_trapz, 1))', ...
                                   sigma(repmat(idxs(j), n_trapz, 1))');
                    sumCdfs = zeros(n_trapz, 1);
                    for k = 1:length(p) - 1
                        sumCdfs = sumCdfs + (1 - cdfs(:, k)) .* prod(cdfs(:, [1:k - 1, k + 1:end]), 2);
                    end
                    y(:, j) = normpdf(x(:, j), means(j), sigma(j)) .* sumCdfs;
                end
                integralsSecondMaximum = trapz(y, 1) .* ((upper_limit - lower_limit) / (n_trapz - 1));

                pHatsOverB = means ./ b;
                
                if ~any(pHatsOverB)
                    revenueWDebias(experiment, i) = 0;
                else                
                    revenueWDebias(experiment, i) = ((integralsMaximum * p') * (integralsSecondMaximum * means')) ...
                                                    / (integralsMaximum * pHatsOverB');
                end
            end
        end

        pTrueTimesB = b .* p;
        maxTrueIdx = datasample(find(pTrueTimesB == max(pTrueTimesB)), 1);
        secondMaxTrueIdx = datasample(find(pTrueTimesB == max(pTrueTimesB([1:maxTrueIdx - 1, maxTrueIdx + 1:end]))), 1);

        figure
        hold

        revenueNoDebiasPlot = 1 - mean(revenueNoDebias) ./ (b(secondMaxTrueIdx) * p(secondMaxTrueIdx));
        revenueDebiasPlot = 1 - mean(revenueDebias) ./ (b(secondMaxTrueIdx) * p(secondMaxTrueIdx));
        revenueWDebiasPlot = 1 - mean(revenueWDebias) ./ (b(secondMaxTrueIdx) * p(secondMaxTrueIdx));

        plot(n_impressions, revenueNoDebiasPlot);
        plot(n_impressions, revenueDebiasPlot);
        plot(n_impressions, revenueWDebiasPlot);

        settingText = num2str(whichSetting);
        save(strcat('sponsoredSearch', settingText, '.mat'), 'n_impressions', 'revenueNoDebias', 'revenueDebias', 'revenueWDebias', ...
                                                                'revenueNoDebiasPlot', 'revenueDebiasPlot', 'revenueWDebiasPlot');

    elseif whichSetting == 3
        p = [0.2 0.15 ones(1, 40) * 0.01]; b = [0.8 1 ones(1, 40)]; n_impressions = [0.5 1 3 5 10 15 20] * 1e3;
        eps = 0.1;
        
        idxs = repmat(1:length(p), length(p), 1);
        idxs = (1 - eye(length(p))) .* idxs';
        idxs(idxs == 0) = [];
        idxs = reshape(idxs, length(p) - 1, length(p));
        idxs = idxs';

        revenueUniform = zeros(n_experiments, length(n_impressions));
        revenueAdaptive = zeros(n_experiments, length(n_impressions));
        revenueWAdaptive = zeros(n_experiments, length(n_impressions));
        for experiment = 1:n_experiments
            fprintf('Experiment: %d\n', experiment);

            for i = 1:length(n_impressions)
                % Biased CTR Learning Phase
                [pHat, means, sigma] = biasedCRT(p, b, n_impressions(i), eps, false);

                % SP Auction Phase
                pHatTimesB = b .* pHat;
                maxHatIdx = datasample(find(pHatTimesB == max(pHatTimesB)), 1);
                secondMaxHatIdx = datasample(find(pHatTimesB == max(pHatTimesB([1:maxHatIdx - 1, maxHatIdx + 1:end]))), 1);

                if pHat(maxHatIdx) == 0
                    revenueUniform(experiment, i) = 0;
                else
                    revenueUniform(experiment, i) = p(maxHatIdx) * b(secondMaxHatIdx) * pHat(secondMaxHatIdx) / pHat(maxHatIdx);
                end

                % Debiased CTR Learning Phase
                [pHatDouble, ~, ~, ~] = debiasedCRT(p, b, n_impressions(i), gamma, true);

                % SP Auction Phase
                pHatTimesB = b .* pHatDouble;
                maxHatIdx = datasample(find(pHatTimesB == max(pHatTimesB)), 1);
                secondMaxHatIdx = datasample(find(pHatTimesB == max(pHatTimesB([1:maxHatIdx - 1, maxHatIdx + 1:end]))), 1);

                if pHatDouble(maxHatIdx) == 0
                    revenueAdaptive(experiment, i) = 0;
                else
                    revenueAdaptive(experiment, i) = p(maxHatIdx) * b(secondMaxHatIdx) ...
                                                       * pHatDouble(secondMaxHatIdx) / pHatDouble(maxHatIdx);
                end

                % Weighted Debiased CTR Learning Phase
                % Same as no selection debiasing!

                % SP Auction Phase
                lower_limit = means - 8 * sigma;
                upper_limit = means + 8 * sigma;
                n_trapz = 1e2;
                x = zeros(n_trapz, length(p));
                y = zeros(size(x));
                for j = 1:length(p)
                    x(:, j) = linspace(lower_limit(j), upper_limit(j), n_trapz);
                    y(:, j) = normpdf(x(:, j), means(j), sigma(j)) .* ...
                            prod(normcdf(repmat(x(:, j), 1, length(p) - 1), ...
                                         means(repmat(idxs(j, :), n_trapz, 1)), ...
                                         sigma(repmat(idxs(j, :), n_trapz, 1))), 2);
                end
                integralsMaximum = trapz(y, 1) .* ((upper_limit - lower_limit) / (n_trapz - 1));
                
                x = zeros(n_trapz, length(p));
                y = zeros(size(x));
                for j = 1:length(p)
                    x(:, j) = linspace(lower_limit(j), upper_limit(j), n_trapz);
                    cdfs = normcdf(repmat(x(:, j), 1, length(p) - 1), ...
                                   means(repmat(idxs(j, :), n_trapz, 1)), ...
                                   sigma(repmat(idxs(j, :), n_trapz, 1)));
                    sumCdfs = zeros(n_trapz, 1);
                    for k = 1:length(p) - 1
                        sumCdfs = sumCdfs + (1 - cdfs(:, k)) .* prod(cdfs(:, [1:k - 1, k + 1:end]), 2);
                    end
                    y(:, j) = normpdf(x(:, j), means(j), sigma(j)) .* sumCdfs;
                end
                integralsSecondMaximum = trapz(y, 1) .* ((upper_limit - lower_limit) / (n_trapz - 1));

                pHatsOverB = means ./ b;
                
                if ~any(pHatsOverB)
                    revenueWAdaptive(experiment, i) = 0;
                else
                    revenueWAdaptive(experiment, i) = ((integralsMaximum * p') * (integralsSecondMaximum * means')) ...
                                                        / (integralsMaximum * pHatsOverB');
                end
            end
        end

        pTrueTimesB = b .* p;
        maxTrueIdx = datasample(find(pTrueTimesB == max(pTrueTimesB)), 1);
        secondMaxTrueIdx = datasample(find(pTrueTimesB == max(pTrueTimesB([1:maxTrueIdx - 1, maxTrueIdx + 1:end]))), 1);

        figure
        hold

        revenueUniformPlot = 1 - mean(revenueUniform) ./ (b(secondMaxTrueIdx) * p(secondMaxTrueIdx));
        revenueAdaptivePlot = 1 - mean(revenueAdaptive) ./ (b(secondMaxTrueIdx) * p(secondMaxTrueIdx));
        revenueWAdaptivePlot = 1 - mean(revenueWAdaptive) ./ (b(secondMaxTrueIdx) * p(secondMaxTrueIdx));

        plot(n_impressions, revenueUniformPlot);
        plot(n_impressions, revenueAdaptivePlot);
        plot(n_impressions, revenueWAdaptivePlot);

        settingText = num2str(whichSetting);
        save(strcat('sponsoredSearch', settingText, '.mat'), 'n_impressions', 'revenueUniform', 'revenueAdaptive', 'revenueWAdaptive', ...
                                                            'revenueUniformPlot', 'revenueAdaptivePlot', 'revenueWAdaptivePlot');
    else
        p = [0.15 0.11 0.1 0.05 0.01]; b1 = linspace(0.8, 1.2, 17); impressions = 10e3;

        idxs = repmat(1:length(p), length(p), 1);
        idxs = (1 - eye(length(p))) .* idxs';
        idxs(idxs == 0) = [];
        idxs = reshape(idxs, length(p) - 1, length(p));
        idxs = idxs';

        v = 1;

        utilitiesNoDebias = zeros(n_experiments, length(b1));
        utilitiesDebias = zeros(n_experiments, length(b1));
        utilitiesWDebias = zeros(n_experiments, length(b1));
        for experiment = 1:n_experiments
            fprintf('Experiment: %d\n', experiment);

            for i = 1:length(b1)
                % Debiased CTR Learning Phase
                b = [b1(i) 0.9 1 2 1];
                [pHat, clicksHistory, means, sigma] = ...
                    debiasedCRT(p, b, impressions, gamma, true);

                % Biased SP Auction
                pHatTimesB = b .* pHat;
                maxHatIdx = datasample(find(pHatTimesB == max(pHatTimesB)), 1);
                secondMaxHatIdx = datasample(find(pHatTimesB == max(pHatTimesB([1:maxHatIdx - 1, maxHatIdx + 1:end]))), 1);

                if pHat(maxHatIdx) == 0
                    utilitiesNoDebias(experiment, i) = 0;
                else
                    utilitiesNoDebias(experiment, i) = v * p(maxHatIdx) - (b(secondMaxHatIdx) * ...
                                                       pHat(secondMaxHatIdx) * p(maxHatIdx)) / pHat(maxHatIdx);
                end

                % Debiased SP Auction
                clickSelection = clicksHistory(1:floor(size(clicksHistory, 1) / 2), :, :);
                clickEstimation = clicksHistory(floor(size(clicksHistory, 1) / 2) + 1:end, :, :);
                pHatSelection = sum(clickSelection(:, :, 1), 1) ./ sum(clickSelection(:, :, 2), 1);
                pHatEstimation = sum(clickEstimation(:, :, 1), 1) ./ sum(clickEstimation(:, :, 2), 1);
                pHatTimesB = b .* pHatSelection;
                maxHatIdx = datasample(find(pHatTimesB == max(pHatTimesB)), 1);
                secondMaxHatIdx = datasample(find(pHatTimesB == max(pHatTimesB([1:maxHatIdx - 1, maxHatIdx + 1:end]))), 1);
                
                if pHatEstimation(maxHatIdx) == 0
                    utilitiesDebias(experiment, i) = 0;
                else                
                    utilitiesDebias(experiment, i) = v * p(maxHatIdx) - (b(secondMaxHatIdx) * ...
                                                     pHatEstimation(secondMaxHatIdx) * p(maxHatIdx)) / pHatEstimation(maxHatIdx);
                end

                % Weighted Debiasing SP Auction
                lower_limit = means - 8 * sigma;
                upper_limit = means + 8 * sigma;
                n_trapz = 1e2;
                x = zeros(n_trapz, length(p));
                y = zeros(size(x));
                for j = 1:length(p)
                    x(:, j) = linspace(lower_limit(j), upper_limit(j), n_trapz);
                    y(:, j) = normpdf(x(:, j), means(j), sigma(j)) .* ...
                            prod(normcdf(repmat(x(:, j), 1, length(p) - 1), ...
                                         means(repmat(idxs(j, :), n_trapz, 1)), ...
                                         sigma(repmat(idxs(j, :), n_trapz, 1))), 2);
                end
                integralsMaximum = trapz(y, 1) .* ((upper_limit - lower_limit) / (n_trapz - 1));

                x = zeros(n_trapz, length(p));
                y = zeros(size(x));
                for j = 1:length(p)
                    x(:, j) = linspace(lower_limit(j), upper_limit(j), n_trapz);
                    cdfs = normcdf(repmat(x(:, j), 1, length(p) - 1), ...
                                   means(repmat(idxs(j, :), n_trapz, 1)), ...
                                   sigma(repmat(idxs(j, :), n_trapz, 1)));
                    sumCdfs = zeros(n_trapz, 1);
                    for k = 1:length(p) - 1
                        sumCdfs = sumCdfs + (1 - cdfs(:, k)) .* prod(cdfs(:, [1:k - 1, k + 1:end]), 2);
                    end
                    y(:, j) = normpdf(x(:, j), means(j), sigma(j)) .* sumCdfs;
                end
                integralsSecondMaximum = trapz(y, 1) .* ((upper_limit - lower_limit) / (n_trapz - 1));

                pHatsOverB = means ./ b;
                
                if ~any(pHatsOverB)
                    utilitiesWDebias(experiment, i) = 0;
                else
                    utilitiesWDebias(experiment, i) = v * (integralsMaximum * p') - ...
                                                      ((integralsSecondMaximum * means') * (integralsMaximum * p')) ...
                                                      / (integralsMaximum * pHatsOverB');
                end
            end
        end

        trueValueIdx = datasample(find(b1 == v), 1);

        figure
        hold

        utilityTrueValue = mean(utilitiesNoDebias(:, trueValueIdx));
        utilitiesNoDebiasPlot = mean(utilitiesNoDebias) ./ utilityTrueValue - 1;
        plot(b1, utilitiesNoDebiasPlot);

        utilityTrueValue = mean(utilitiesDebias(:, trueValueIdx));
        utilitiesDebiasPlot = mean(utilitiesDebias) ./ utilityTrueValue - 1;
        plot(b1, utilitiesDebiasPlot);

        utilityTrueValue = mean(utilitiesWDebias(:, trueValueIdx));
        utilitiesWDebiasPlot = mean(utilitiesWDebias) ./ utilityTrueValue - 1;
        plot(b1, utilitiesWDebiasPlot);

        settingText = num2str(whichSetting);
        save(strcat('sponsoredSearch', settingText, '.mat'), 'b1', 'utilitiesNoDebias', 'utilitiesDebias', 'utilitiesWDebias', ...
                                                'utilitiesNoDebiasPlot', 'utilitiesDebiasPlot', 'utilitiesWDebiasPlot');
    end
end
