function [clicks, means, sigma] = crtLearning(n_actions, n_trials, p)

n_updates = 1;
current_action = 1;
clicks = zeros(n_trials, n_actions);

for i = 1:n_trials * n_actions
    clicks(n_updates, current_action) = computeReward(p(current_action));

    current_action = current_action + 1;
    if(current_action == n_actions + 1)
        current_action = 1;
        n_updates = n_updates + 1;
    end
end

means = sum(clicks) / n_updates;
variance = (means - means.^2) / n_updates; % This is possible because clicks = 0 or 1
variance(variance < 0) = 0;
sigma = sqrt(variance);
sigma(sigma < 1e-4) = 1e-4;

end
