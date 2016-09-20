function [max_Q, r] = QLearning(n_states, n_actions, n_experiments, start, goal, n_cols, n_steps, gamma, exp, reward_array, policy)
% Q-Learning algorithm

max_Q = zeros(n_experiments, n_steps);
r = zeros(n_experiments, n_steps);

for experiment = 1:n_experiments
    Q = zeros(n_states, n_actions);
    current_state = start;
    n_alpha = zeros(n_states, n_actions);
    n_eps = zeros(n_states);
    
    fprintf('Experiment: %d\n', experiment);
    
    for i = 1:n_steps
        n_eps(current_state) = n_eps(current_state) + 1;
        
        [action, reward, next_state] = step(current_state, Q, n_states, n_cols, n_eps, start, goal, experiment, i, reward_array, policy);
        
        n_alpha(current_state, action) = n_alpha(current_state, action) + 1;

        if(current_state ~= goal)
            delta = reward + gamma * max(Q(next_state, :)) - Q(current_state, action);
        else
            delta = reward - Q(current_state, action);
        end
        Q(current_state, action) = Q(current_state, action) + (1 / n_alpha(current_state, action)^exp) * delta;

        if(mod(i, 5000) == 0)
            fprintf('Trial: %d\n', i);
        end
        
        r(experiment, i) = reward;
        max_Q(experiment, i) = max(Q(start, :));
        current_state = next_state;
    end
end

end

