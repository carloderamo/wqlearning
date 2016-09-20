function [max_Q, r] = WeightedQLearning(n_states, n_actions, n_experiments, start, goal, n_cols, n_steps, gamma, exp, reward_array, policy)
% Weighted Q-Learning

max_Q = zeros(n_experiments, n_steps);
r = zeros(n_experiments, n_steps);

idxs = repmat(1:n_actions, n_actions, 1);
idxs = (1 - eye(n_actions)) .* idxs';
idxs(idxs == 0) = [];
idxs = reshape(idxs, n_actions - 1, n_actions);
idxs = idxs';

for experiment = 1:n_experiments
    Q = zeros(n_states, n_actions);
    Q2 = Q;
    weightsVar = Q;
    current_state = start;
    n_alpha = zeros(n_states, n_actions);
    n_eps = zeros(n_states);
    sigma = ones(size(Q)) * 1e10;
    probs = ones(1, n_actions) * (1 / n_actions);
	n_updates = zeros(n_states, n_actions);
    
    fprintf('Experiment: %d\n', experiment);
    
    for i = 1:n_steps
        n_eps(current_state) = n_eps(current_state) + 1;
        
        [action, reward, next_state] = step(current_state, Q, n_states, n_cols, n_eps, start, goal, experiment, i, reward_array, policy, probs);

        means = Q(next_state, :);

        current_sigma = sigma(next_state, :);
	
        current_sigma(current_sigma < 1e-4) = 1e-4;
        
%         lower_limit = means - 8 * current_sigma;
%         upper_limit = means + 8 * current_sigma;
%         n_trapz = 1e2;
%         x = zeros(n_trapz, n_actions);
%         y = zeros(size(x));
%         for j = 1:n_actions
%            x(:, j) = linspace(lower_limit(j), upper_limit(j), n_trapz);
%            y(:, j) = normpdf(x(:, j), means(j), current_sigma(j)) .* ...
%                    prod(normcdf(repmat(x(:, j), 1, n_actions - 1), ...
%                                 means(repmat(idxs(j, :), n_trapz, 1)), ...
%                                 current_sigma(repmat(idxs(j, :), n_trapz, 1))), 2);
%         end
%         probs = trapz(y, 1) .* ((upper_limit - lower_limit) / (n_trapz - 1));
%         
%         fprintf('\nIntegral: ');
%         sum(probs)
%         W = probs * Q(next_state, :)';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% SAMPLING PROBS %%%%%%
        n_samples = 1000;
        samples = repmat(means, n_samples, 1) + repmat(current_sigma, n_samples, 1) .* randn(n_samples, n_actions);
        [~, max_idxs] = max(samples');

        occCount = zeros(1, 4);
        [occ, val] = hist(max_idxs, unique(max_idxs));
        occCount(val(occ > 0)) = occ(occ > 0);
        probs = occCount / n_samples;

        W = probs * Q(next_state, :)';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% PROBABILITIES PRODUCT %%%%%
%         F = zeros(n_actions, n_actions - 1);
%         for j = 1:n_actions
%             diffMeans = means(idxs(j, :)) - means(j);
%             diffSigma = sqrt(current_sigma(idxs(j, :)).^2 + current_sigma(j)^2);
%             F(j, :) = normcdf(0, diffMeans, diffSigma);
%         end
%         prods = prod(F, 2);
%         
%         fprintf('Prod: ');
%         sum(prods)
%         W = prods' * Q(next_state, :)'
%         fprintf('...');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if(current_state ~= goal)
            target = reward + gamma * W;
        else
            target = reward;
        end
        
        n_alpha(current_state, action) = n_alpha(current_state, action) + 1;

		n_updates(current_state, action) = n_updates(current_state, action) + 1;
        
        alpha = 1 / n_alpha(current_state, action)^exp;
        Q(current_state, action) = (1 - alpha) * Q(current_state, action) + alpha * target;
        Q2(current_state, action) = (1 - alpha) * Q2(current_state, action) + alpha * target^2;
        
		if n_updates(current_state, action) > 1
            weightsVar(current_state, action) = (1 - alpha)^2 * weightsVar(current_state, action) + alpha^2;
		    n = 1 / weightsVar(current_state, action);
		    diff = Q2(current_state, action) - Q(current_state, action)^2;
		    if diff < 0
		        diff = 0;
		    end
		    sigma(current_state, action) = sqrt(diff / n);
		end

        if(mod(i, 500) == 0)
            fprintf('Trial: %d\n', i);
        end

        r(experiment, i) = reward;
        max_Q(experiment, i) = max(Q(start, :));
        current_state = next_state;
    end
end

end
