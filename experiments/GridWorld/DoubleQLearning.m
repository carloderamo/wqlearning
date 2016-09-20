function [max_Q, r] = DoubleQLearning(n_states, n_actions, n_experiments, start, goal, n_cols, n_steps, gamma, exp, reward_array, policy)
% Double Q-Learning algorithm

max_Q = zeros(n_experiments, n_steps);
r = zeros(n_experiments, n_steps);

which_Q = uint8(rand(n_experiments, n_steps)) + 1;

for experiment = 1:n_experiments
    Q = zeros(n_states, n_actions, 2);
    current_state = start;
    n_alpha = zeros(n_states, n_actions, 2);
    n_eps = zeros(n_states);
    
    fprintf('Experiment: %d\n', experiment);
    
    for i = 1:n_steps
        n_eps(current_state) = n_eps(current_state) + 1;
        
        [action, reward, next_state] = step(current_state, Q, n_states, n_cols, n_eps, start, goal, experiment, i, reward_array, policy);
        
        n_alpha(current_state, action, which_Q(experiment, i)) = n_alpha(current_state, action, which_Q(experiment, i)) + 1;
        
        [~, argmax] = max(Q(next_state, :, which_Q(experiment, i)));

        if(current_state ~= goal)
            delta = reward + gamma * Q(next_state, argmax, 3 - which_Q(experiment, i)) - Q(current_state, action, which_Q(experiment, i));
        else
            delta = reward - Q(current_state, action, which_Q(experiment, i));
        end
        Q(current_state, action, which_Q(experiment, i)) = Q(current_state, action, which_Q(experiment, i)) ...
            + (1 / n_alpha(current_state, action, which_Q(experiment, i))^exp) * delta;

        if(mod(i, 5000) == 0)
            fprintf('Trial: %d\n', i);
        end
        
        r(experiment, i) = reward;
        mean_Q = mean(Q, 3);
        max_Q(experiment, i) = max(mean_Q(start, :));
        current_state = next_state;
    end
end

end
