clf
hold

x = 0.02:0.01:0.1;

mseBar = [reshape(bias'.^2, numel(bias), 1) reshape(variance', numel(bias), 1)];

set(gca,'FontSize', 22)
pbaspect([1 1 1])

width = 0.25;
gap = 0.002;

me = mseBar(1:4:size(mseBar, 1) - 3, :);
h = bar(x - gap, me, width, 'stacked');

xlabel('Max probability in range');


set(h(1), 'FaceColor', [0 0 0.4])
set(h(2), 'FaceColor', [0.5 0.5 1])

de =  mseBar(2:4:size(mseBar, 1) - 2, :);
h = bar(x, de, width, 'stacked');

xlabel('Max probability in range');

ylim([0, 2e-5])
xlim([0.01,0.11])

set(h(1), 'FaceColor', [0.4 0 0])
set(h(2), 'FaceColor', [1 0 0])

mme = mseBar(3:4:size(mseBar, 1) - 1, :);
h = bar(x + gap, mme, width, 'stacked');

xlabel('Max probability in range');


set(h(1), 'FaceColor', [1 0.65 0])
set(h(2), 'FaceColor', [1 1 0])

h = bar(x + 2 * gap, mseBar(4:4:size(mseBar, 1), :), width, 'stacked');

xlabel('Max probability in range');

set(h(1), 'FaceColor', [0 0.4 0])
set(h(2), 'FaceColor', [0 1 0])

saveas(gca,'internetAds-probs.pdf');
