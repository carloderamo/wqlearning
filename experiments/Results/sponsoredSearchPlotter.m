subplot(2, 2, 1)
hold

plot(f1.n_impressions, f1.revenueNoDebiasPlot, 'b')
plot(f1.n_impressions, f1.revenueDebiasPlot, 'r')
plot(f1.n_impressions, f1.revenueWDebiasPlot, 'g')

xlabel('Rounds of Exploration')
ylabel('Expected Revenue Loss')
axis([0 15000 -0.02 0.1])
legend('No Selection Debiasing', 'Double Selection Debiasing', 'WG Selection Debiasing')

subplot(2, 2, 2)
hold

plot(f2.n_impressions, f2.revenueNoDebiasPlot, 'b')
plot(f2.n_impressions, f2.revenueDebiasPlot, 'r')
plot(f2.n_impressions, f2.revenueWDebiasPlot, 'g')

xlabel('Rounds of Exploration')
ylabel('Expected Revenue Loss')
axis([0 15000 -0.02 0.05])
legend('No Selection Debiasing', 'Double Selection Debiasing', 'WG Selection Debiasing')

subplot(2, 2, 3)
hold

plot(f3.n_impressions, f3.revenueUniformPlot, 'b')
plot(f3.n_impressions, f3.revenueAdaptivePlot, 'r')
plot(f3.n_impressions, f3.revenueWAdaptivePlot, 'g')

xlabel('Rounds of Exploration')
ylabel('Expected Revenue Loss')
axis([0 20000 0 0.6])
legend('Uniform', 'Adaptive with Debiasing', 'WG Adaptive with Debiasing')

subplot(2, 2, 4)
hold

plot(f4.b1, f4.utilitiesNoDebiasPlot, 'b')
plot(f4.b1, f4.utilitiesDebiasPlot, 'r')
plot(f4.b1, f4.utilitiesWDebiasPlot, 'g')

xlabel('Bid Price')
ylabel('Player Utility Gain')
axis([0.8 1.2 -0.2 0.1])
legend('No ELM Debiasing', 'Double ELM Debiasing', 'WG ELM Debiasing')

x = [1 1];
y = [-5 5];
plot(x, y, '--')
