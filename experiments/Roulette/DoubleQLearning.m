function [avg_Q] = DoubleQLearning(n_actions, n_experiments, n_trials, gamma, exp)
% Double Q-Learning algorithm

avg_Q = zeros(n_experiments, n_trials);

parfor experiment = 1:n_experiments
    Q = zeros(2, n_actions);
    prevQ = Q;
    n_updates = 0;
    current_action = 1;
    n_lambda = ones(2, n_actions);
    
    which_Q = int8(rand);
    other_Q = 2 - which_Q;
    which_Q = which_Q + 1;
    
    [~, argmax] = max(prevQ(which_Q, :));
    
    fprintf('Experiment: %d\n', experiment);
    
    for i = 1:n_trials * n_actions
        reward = computeReward(current_action);

        delta = reward + gamma * prevQ(other_Q, argmax) - Q(which_Q, current_action);
        Q(which_Q, current_action) = Q(which_Q, current_action) + (1 / n_lambda(which_Q, current_action)^exp) * delta;

        n_lambda(which_Q, current_action) = n_lambda(which_Q, current_action) + 1;

        current_action = current_action + 1;
        if(current_action == n_actions + 1)
            current_action = 1;
            which_Q = int8(rand);
            other_Q = 2 - which_Q;
            which_Q = which_Q + 1;
            
            [~, argmax] = max(prevQ(which_Q, :));
            
            prevQ = Q;
            
            n_updates = n_updates + 1;
            avg_Q(experiment, n_updates) = mean(mean(prevQ));
            
            if(mod(n_updates, 5000) == 0)
                fprintf('Trial: %d\n', n_updates);
            end            
        end
    end
end

end
