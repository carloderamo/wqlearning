%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script aims to replicate the grid world     %
% experiments presented in H. V. Hasselt's article %
% "Double Q-Learning".                             %
%                                                  %
% Written by: Carlo D'Eramo                        %
%                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

gamma = 0.95;
exp = input('Exp: ');
n_experiments = 10000;
n_actions = 4;
n_steps = 10000;
reward_type = input('Reward Type: '); % use 'g' for gaussian reward or 'ng' for reward described by the paper
policy = input('Policy: '); % use 'neps' for W policy or 'eps for epsilon greedy

n_rows = 3;
n_cols = 3;
start = 7;
goal = 3;
n_states = n_rows * n_cols;

r1 = 0; r2 = 0; r3 = 0;
max_Q1 = 0; max_Q2 = 0; max_Q3 = 0;

% if strcmp(policy, 'eps') || strcmp(policy, 'epsAndGreedy')
% 	if(strcmp(reward_type, 'g'))
% 		reward_array = normrnd(ones(n_experiments, n_steps) * -1, ones(n_experiments, n_steps) * 1);
%     elseif(strcmp(reward_type, 'ng'))
% 		reward_array = rand(n_experiments, n_steps);
% 		reward_array(reward_array < 0.5) = -12;
% 		reward_array(reward_array >= 0.5) = 10;
% 	end
% 
% 	figure
% 	[max_Q1, r1] = QLearning(n_states, n_actions, n_experiments, start, goal, n_cols, n_steps, gamma, exp, reward_array, policy);
% 	plot(mean(max_Q1));
% 	figure
% 	plot(mean(r1));
% 
% 	if(strcmp(reward_type, 'g'))
% 		reward_array = normrnd(ones(n_experiments, n_steps) * -1, ones(n_experiments, n_steps) * 1);
%     elseif(strcmp(reward_type, 'ng'))
% 		reward_array = rand(n_experiments, n_steps);
% 		reward_array(reward_array < 0.5) = -12;
% 		reward_array(reward_array >= 0.5) = 10;
% 	end
% 	figure
% 	[max_Q2, r2] = DoubleQLearning(n_states, n_actions, n_experiments, start, goal, n_cols, n_steps, gamma, exp, reward_array, policy);
% 	plot(mean(max_Q2));
% 	figure
% 	plot(mean(r2));
% end
if(strcmp(reward_type, 'g'))
	reward_array = normrnd(ones(n_experiments, n_steps) * -1, ones(n_experiments, n_steps) * 1);
elseif(strcmp(reward_type, 'ng'))
	reward_array = rand(n_experiments, n_steps);
	reward_array(reward_array < 0.5) = -12;
	reward_array(reward_array >= 0.5) = 10;
end

figure
[max_Q3, r3] = W(n_states, n_actions, n_experiments, start, goal, n_cols, n_steps, gamma, exp, reward_array, policy);
plot(mean(max_Q3));
figure
plot(mean(r3));

exp_text = num2str(exp);
reward_type_text = strcat('-', reward_type);
policy_text = strcat('-', policy);
save(strcat('grid-exp', exp_text, reward_type_text, policy_text, '.mat'), 'r1', 'r2', 'r3', 'max_Q1', 'max_Q2', 'max_Q3');
