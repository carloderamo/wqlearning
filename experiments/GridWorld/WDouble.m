function [max_Q, r] = WDouble(n_states, n_actions, n_experiments, start, goal, n_cols, n_steps, gamma, exp, reward_array, policy)
% Weighted Double Q-Learning

max_Q = zeros(n_experiments, n_steps);
r = zeros(n_experiments, n_steps);

which_Q = uint8(rand(n_experiments, n_steps)) + 1;

for experiment = 1:n_experiments
    Q = zeros(n_states, n_actions, 2);
    Q2 = Q;
    weightsVar = Q;
    current_state = start;
    n_alpha = zeros(n_states, n_actions, 2);
    n_eps = zeros(n_states);
    sigma = ones(size(Q)) * 1e10;
    probs = ones(1, n_actions, 2) * (1 / n_actions);
	n_updates = zeros(n_states, n_actions, 2);
    
    fprintf('Experiment: %d\n', experiment);
    
    for i = 1:n_steps
        n_eps(current_state) = n_eps(current_state) + 1;
        
        [action, reward, next_state] = step(current_state, Q, n_states, ...
                                            n_cols, n_eps, start, goal, ...
                                            experiment, i, reward_array, ...
                                            policy, probs(:, :, which_Q(experiment, i)));
        
        means = Q(next_state, :, which_Q(experiment, i));
        current_sigma = sigma(next_state, :, which_Q(experiment, i));
        current_sigma(current_sigma < 1e-4) = 1e-4;
        
        n_samples = 1000;
        samples = repmat(means, n_samples, 1) + repmat(current_sigma, n_samples, 1) .* randn(n_samples, n_actions);
        [~, max_idxs] = max(samples');

        occCount = zeros(1, 4);
        [occ, val] = hist(max_idxs, unique(max_idxs));
        occCount(val(occ > 0)) = occ(occ > 0);
        probs(:, :, which_Q(experiment, i)) = occCount / n_samples;

        W = probs(:, :, which_Q(experiment, i)) * Q(next_state, :, 3 - which_Q(experiment, i))';
        
        if(current_state ~= goal)
            target = reward + gamma * W;
        else
            target = reward;
        end
        
        n_alpha(current_state, action, which_Q(experiment, i)) = n_alpha(current_state, action, which_Q(experiment, i)) + 1;
        
		n_updates(current_state, action, which_Q(experiment, i)) = n_updates(current_state, action, which_Q(experiment, i)) + 1;
        
        alpha = 1 / n_alpha(current_state, action, which_Q(experiment, i))^exp;
        Q(current_state, action, which_Q(experiment, i)) = ...
            (1 - alpha) * Q(current_state, action, which_Q(experiment, i)) + alpha * target;
        Q2(current_state, action, which_Q(experiment, i)) = ...
            (1 - alpha) * Q2(current_state, action, which_Q(experiment, i)) + alpha * target^2;
        
		if n_updates(current_state, action) > 1
            weightsVar(current_state, action, which_Q(experiment, i)) = ...
                (1 - alpha)^2 * weightsVar(current_state, action, which_Q(experiment, i)) + alpha^2;
		    n = 1 / weightsVar(current_state, action, which_Q(experiment, i));
		    diff = Q2(current_state, action, which_Q(experiment, i)) - Q(current_state, action, which_Q(experiment, i))^2;
		    if diff < 0
		        diff = 0;
		    end
		    sigma(current_state, action, which_Q(experiment, i)) = sqrt(diff / n);
		end

        if(mod(i, 500) == 0)
            fprintf('Trial: %d\n', i);
        end

        r(experiment, i) = reward;
        max_Q(experiment, i) = max(mean(Q(start, :, :), 3));
        current_state = next_state;
    end

end