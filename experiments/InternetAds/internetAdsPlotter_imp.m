clf
hold

x = 1000:1000:10000;
mseBar = [reshape(bias'.^2, numel(bias), 1) reshape(variance', numel(bias), 1)];

set(gca,'FontSize', 22)
pbaspect([1 1 1])

width = 0.25;
gap = 180;

me = mseBar(1:4:size(mseBar, 1) - 3, :);
h = bar(x - gap, me, width, 'stacked');

xlabel('N. impressions');
ylabel('MSE')

legend('ME bias²', 'ME variance')

set(h(1), 'FaceColor', [0 0 0.4])
set(h(2), 'FaceColor', [0.5 0.5 1])

de =  mseBar(2:4:size(mseBar, 1) - 2, :);
h = bar(x, de, width, 'stacked');

xlabel('N. impressions');
ylabel('MSE')
ylim([0, 8e-5])

set(h(1), 'FaceColor', [0.4 0 0])
set(h(2), 'FaceColor', [1 0 0])

mme = mseBar(3:4:size(mseBar, 1) - 1, :);
h = bar(x + gap, mme, width, 'stacked');

xlabel('N. impressions');
ylabel('MSE')

set(h(1), 'FaceColor', [1 0.65 0])
set(h(2), 'FaceColor', [1 1 0])

h = bar(x + 2*gap, mseBar(4:4:size(mseBar, 1), :), width, 'stacked');

xlabel('N. impressions');
ylabel('MSE')

set(h(1), 'FaceColor', [0 0.4 0])
set(h(2), 'FaceColor', [0 1 0])

legend('ME bias²', 'ME variance', 'DE bias²', 'DE variance','MME bias²', 'MME variance', 'WE bias²', 'WE variance');

saveas(gca,'internetAds-impressions.pdf');
