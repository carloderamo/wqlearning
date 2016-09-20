function [max_Q, r] = BiasCorrectedQLearning(n_states, n_actions, n_experiments, start, goal, n_cols, n_steps, gamma, exp, reward_array, policy)
% BiasCorrectedQ-Learning algorithm

max_Q = zeros(n_experiments, n_steps);
r = zeros(n_experiments, n_steps);
euler=0.5772;

for experiment = 1:n_experiments
    Q = zeros(n_states, n_actions);
    current_state = start;
    nVisits = ones(n_states, n_actions);
    B = zeros(n_states, n_actions);
	Rvar=zeros(n_states, n_actions);      %%% corretta inizializzazione a 0 ??????? ( nella prima iterazione la Q non viene corretta da B)
	Rmean=zeros(n_states,n_actions);    
	n_eps = zeros(n_states);
    
    fprintf('Experiment: %d\n', experiment);
    
    for i = 1:n_steps
        n_eps(current_state) = n_eps(current_state) + 1;
        
        [action, reward, next_state] = stepForBiasCQL(current_state, Q, n_states, n_cols, n_eps, start, goal, experiment, i, reward_array, policy);
        
        prevMean=Rmean(current_state, action);
		prevVar=Rvar(current_state, action);
		prevSigma=sqrt(prevVar/nVisits(current_state,action));

		Rmean(current_state, action) = prevMean + (reward - prevMean)/nVisits(current_state, action);
		Rvar(current_state, action) = (prevVar + (reward- prevMean)*(reward - Rmean(current_state,action)))/nVisits(current_state,action);
		
		bM=(2*log(n_actions) - log(log(n_actions)) - log(4*pi))^0.5;
		B(current_state, action)=(euler/bM + bM)*prevSigma;
		
        if(current_state ~= goal)
            delta = reward + gamma * max(Q(next_state, :)) - Q(current_state, action);
        else
            delta = reward - Q(current_state, action);
        end
        Q(current_state, action) = Q(current_state, action) + (1 / nVisits(current_state, action)^exp) * (delta - B(current_state, action));

        if(mod(i, 5000) == 0)
            fprintf('Trial: %d\n', i);
        end
        		nVisits(current_state, action) = nVisits(current_state, action) + 1;

        r(experiment, i) = reward;
        max_Q(experiment, i) = max(Q(start, :));
        current_state = next_state;
    end
end

end

