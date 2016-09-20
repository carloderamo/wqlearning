%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script aims to replicate the roulette       %
% experiments presented in H. V. Hasselt's article %
% "Double Q-Learning".                             %
%                                                  %
% Written by: Carlo D'Eramo                        %
%                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

gamma = 0.95;
exp = input('Exp: ');
n_experiments = 10;
n_actions = 158;
n_trials = 1e5;

avg_Q1 = 0;
avg_Q2 = 0;
avg_Q3 = 0;

avg_Q1 = QLearning(n_actions, n_experiments, n_trials, gamma, exp);
plot(mean(avg_Q1));

figure
avg_Q2 = DoubleQLearning(n_actions, n_experiments, n_trials, gamma, exp);
plot(mean(avg_Q2));

figure
avg_Q3 = W(n_actions, n_experiments, n_trials, gamma, exp);
plot(mean(avg_Q3));

exp_text = num2str(exp);
save(strcat('roulette-exp', exp_text,'.mat'), 'avg_Q1', 'avg_Q2', 'avg_Q3');
