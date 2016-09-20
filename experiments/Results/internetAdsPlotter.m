clf
hold

x = 0.02:0.01:0.1;
mseBar = [reshape(f.bias'.^2, numel(f.bias), 1) reshape(f.variance', numel(f.bias), 1)];

h = bar(x - 0.002, mseBar(1:3:size(mseBar, 1) - 2, :), 0.2, 'stacked');

xlabel('Max probability in range');
ylabel('MSE')

legend('ME bias²', 'ME variance')

set(h(1), 'FaceColor', [0 0 0.4])
set(h(2), 'FaceColor', [0.5 0.5 1])

h = bar(x, mseBar(2:3:size(mseBar, 1) - 1, :), 0.2, 'stacked');

xlabel('Max probability in range');
ylabel('MSE')

set(h(1), 'FaceColor', [0.4 0 0])
set(h(2), 'FaceColor', [1 0 0])

h = bar(x + 0.002, mseBar(3:3:size(mseBar, 1), :), 0.2, 'stacked');

xlabel('Max probability in range');
ylabel('MSE')

set(h(1), 'FaceColor', [0 0.4 0])
set(h(2), 'FaceColor', [0 1 0])

legend('ME bias²', 'ME variance', 'CV bias²', 'CV variance', 'WE bias²', 'WE variance')
