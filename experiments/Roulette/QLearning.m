function [avg_Q] = QLearning(n_actions, n_experiments, n_trials, gamma, exp)
% Q-Learning algorithm

avg_Q = zeros(n_experiments, n_trials);

parfor experiment = 1:n_experiments
    Q = zeros(1, n_actions);
    prevQ = Q;
    n_updates = 0;
    current_action = 1;
    n_lambda = ones(n_actions);
    
    fprintf('Experiment: %d\n', experiment);
    
    for i = 1:n_trials * n_actions
        reward = computeReward(current_action);
        delta = reward + gamma * max(prevQ) - Q(current_action);
        Q(current_action) = Q(current_action) + (1 / n_lambda(current_action)^exp) * delta;
        
        n_lambda(current_action) = n_lambda(current_action) + 1;

        current_action = current_action + 1;
        if(current_action == n_actions + 1)
            current_action = 1;
            prevQ = Q;
            
            n_updates = n_updates + 1;
            avg_Q(experiment, n_updates) = mean(prevQ);
            
            if(mod(n_updates, 5000) == 0)
                fprintf('Trial: %d\n', n_updates);
            end
        end
    end
end

end

