subplot(2,1,1)
hold

plot(mean(f1.avg_Q1), 'b')
plot(mean(f1.avg_Q2), 'r')
plot(mean(f1.avg_Q3), 'g')

xlabel('Number of Trials')
ylabel('Expected Profit')
axis([0 500000 -1 41])
legend('Q-Learning', 'Double Q-Learning', 'WG Q-Learning')

x = [0 500000];
y = [-0.053 -0.053];
plot(x, y, '--')

subplot(2,1,2)
hold

plot(mean(f2.avg_Q1), 'b')
plot(mean(f2.avg_Q2), 'r')
plot(mean(f2.avg_Q3), 'g')

xlabel('Number of Trials')
ylabel('Expected Profit')
axis([0 500000 -1 41])
legend('Q-Learning', 'Double Q-Learning', 'WG Q-Learning')

x = [0 500000];
y = [-0.053 -0.053];
plot(x, y, '--')
