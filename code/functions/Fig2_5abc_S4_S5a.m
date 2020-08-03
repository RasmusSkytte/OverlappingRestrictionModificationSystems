close all
clearvars

% Load data
load('../data/Fig2/dataNitish.mat');

% Group by 'genus'
g = 3;
s = S(:, g);
labels = strrep(uniqueLabels{g}, '/', '_');
labels = cellfun(@(s) ['_' s], labels, 'UniformOutput', false);

% Prepare figure folders
if ~exist('../../figures/Figure_2', 'dir')
    mkdir('../../figures/Figure_2')
end

if ~exist('../../figures/Figure_5', 'dir')
    mkdir('../../figures/Figure_5')
end

% Prepare figures
fh1 = figure('Resize', 'off'); clf; 
ax1 = axes; box on;

fh2 = figure('Resize', 'off'); clf; 
ax2 = axes; box on;

fh3 = figure('Resize', 'off'); clf; 
ax3 = axes; box on;

fh4 = figure('Resize', 'off'); clf; 
axNet = axes;

fh5 = figure('Resize', 'off'); clf; 
ax5 = axes; hold on; box on;


fh4.Color = [1 1 1];
set(fh4, 'PaperPositionMode', 'auto')
set(fh4, 'InvertHardcopy', 'off')

% Update figure size
drawnow;
fh1.Position = [488 342 560 420];
fh2.Position = [488 342 560 420];
fh3.Position = [488 342 560 420];
fh5.Position = [488 342 560 420];


fh4.Position(2) = 200;
fh4.Position(4) = fh4.Position(3);

% Ensure the data is sorted
[s, sortIndex] = sort(s);
S = S(sortIndex, :);
A_ij = A_ij(sortIndex, :);


% Allocate arrays
nRM = [];
o_ij = [];
max_o_ij = [];
nSamples = 0;
nUnviable = 0;
nViable = 0;

avgOverlap = [];
avgRM      = [];

% Loop over subsets to generate figures
for d = 1:numel(unique(s))
    
    % Locate the start and size of the subset
    is = find(s == d, 1);
    ns = sum(s == d);
   
    % Extract the subset
    A_ij_s = A_ij(is:is+(ns-1), :);
    
    % Filter the subset for unique combinations
    A_ij_s = unique(A_ij_s, 'rows');
    
    % Count the samples
    nSamples = nSamples + size(A_ij_s, 1);
        
    if size(A_ij_s, 1) < 20
        nUnviable = nUnviable + 1;
        continue
    end
    
    % Count the viable genera
    nViable = nViable + 1;
    
    % Label the adjacency matrix
    A_ij_s = A_ij_s .* repmat(1:size(A_ij_s, 2), size(A_ij_s, 1), 1);
    
    % Convert to list of RM format
    Bs = cell(size(A_ij_s, 1), 1);
    for i = 1:size(A_ij_s, 1)
        Bs{i} = A_ij_s(i, A_ij_s(i, :) > 0);
    end
    

    % Perform sanity checks
    if numel(unique(S(is:is+(ns-1), 2))) > 1
        warning('The genus "%s" has several sup-groups', uniqueLabels{g}{d})
    end
    if numel(unique(S(is:is+(ns-1), 1))) > 1
        warning('The genus "%s" has several groups', uniqueLabels{g}{d})
    end
    
    % Get the distribution of number of RM systems
    nRM_s = sum(A_ij_s > 0, 2);
    nRM = [nRM; nRM_s];
    
    % Compute relative overlap between species
    A_ij_s = A_ij_s > 0;
    o_ij_s = nan(numel(Bs), numel(Bs));
    for i = 1:numel(Bs)
        for j = 1:numel(Bs)
            if i == j
                o_ij_s(i, j) = nan;
                continue;
            end
            o_ij_s(i, j) = sum(and(A_ij_s(i, :), A_ij_s(j, :))) / sum( or(A_ij_s(i, :), A_ij_s(j, :)));
        end
    end
    
    % Store values for histograms
    o_ij = [o_ij; o_ij_s(:)];
    max_o_ij = [max_o_ij; max(o_ij_s)'];
    
    % Plot the RM distribution
    histogram(ax1, nRM, -0.5:(max(nRM)+0.5), 'LineWidth', 1, 'FaceAlpha', 0.4, 'Normalization', 'Probability');
    
    % Plot the histogram of overlap
    histogram(ax2, o_ij(:), 0:0.1:1, 'LineWidth', 1, 'FaceAlpha', 0.4, 'Normalization', 'Probability');
    
    % Plot the histogram of min overlap
    histogram(ax3, max_o_ij, 0:0.1:1, 'LineWidth', 1, 'FaceAlpha', 0.4, 'Normalization', 'Probability');
    

    % Draw the graph at the end of the run
    [~, A_RM_s] = computeNetwork(Bs);

    if isempty(A_RM_s) % No RM systems present
        continue;
    end

    % Compute overlap
    overlap_s = (arrayfun(@(r) sum([Bs{:}]==r), unique([Bs{:}])) - 1)/(size(A_ij_s, 1) - 1);


    % Plot the RM network
    G = graph(A_RM_s);
    l = arrayfun(@(d) num2str(d), find(sum(A_ij_s, 1) > 0)', 'UniformOutput', false);
    l(G.degree < quantile(G.degree, 0.7)) = {' '};
    plotNetwork(axNet, G, l)

    % Adjust axes
    axNet.Position = [0 0 1 1];
    axNet.Visible = 'off';
    axNet.XLim = [-1.4 1.4];
    axNet.YLim = [-1.4 1.4]+0.2;

    % Label the figure
    text(axNet, -1.3,  1.50, labels{d}(2:end), 'FontSize', 18, 'FontWeight', 'Bold');
    text(axNet, -1.3,  1.35, sprintf('# RM = %d,', size(A_RM_s, 1)), 'FontSize', 14);
    text(axNet, -1.3,  1.25, sprintf('# samples = %d', size(A_ij_s, 1)), 'FontSize', 14);
    text(axNet, -0.4,  1.35, sprintf('{\\langle}#RM{\\rangle} = %.2f,', mean(sum(A_ij_s > 0, 2))), 'FontSize', 14);
    text(axNet, -0.4,  1.25, sprintf('{\\langle}p{\\rangle} = %.4f,', mean(overlap_s)), 'FontSize', 14);
    
    % Save the network figure
    print(fh4, sprintf('../../figures/Figure_5/fig5%s.tif', labels{d}), '-dtiff')

    % Report overlap
    fprintf('<p_{%s}> = %.4f, n = %d\n', uniqueLabels{g}{d}, mean(overlap_s), size(A_ij_s, 1))

    % Store average measures
    avgOverlap(end + 1) = mean(overlap_s);
    avgRM(end + 1)      = mean(sum(A_ij_s > 0, 2));
    
end

% Adjust the figures
xlabel(ax1, '# RM')
ylabel(ax1, 'pmf')
ax1.XLim = [-0.5 max(arrayfun(@(ch) ch.BinEdges(end), ax1.Children) + 1)];
ax1.YLim = [0 0.5];
ax1.XTick = 0:2:24;
ax1.Position = [0.165 0.18 0.83 0.79];
ax1.FontSize = 20;
ax1.LineWidth = 1;
ytickformat(ax1, '%.1f');

xlabel(ax2, '# Intersection / # Union')
ylabel(ax2, 'pmf')
ax2.XLim = [-0.1 1.1];
ax2.YLim = [0 1];
ax2.Position = [0.165 0.18 0.83 0.79];
ax2.FontSize = 20;
ax2.LineWidth = 1;
ytickformat(ax2, '%.1f');

xlabel(ax3, 'max(# Intersection / # Union)')
ylabel(ax3, 'pmf')
ax3.XLim = [-0.1 1.1];
ax3.YLim = [0 1];
ax3.Position = [0.165 0.18 0.83 0.79];
ax3.FontSize = 20;
ax3.LineWidth = 1;
ytickformat(ax3, '%.1f');

scatter(ax5, avgOverlap, avgRM, 'filled')
xlabel(ax5, '{\langle}p\rangle')
ylabel(ax5, '{\langle}#RM\rangle')
ax5.Position = [0.165 0.19 0.77 0.78];
ax5.FontSize = 20;
ax5.LineWidth = 1;
ax5.XLim = [0 0.025];
ax5.YLim = [1 4];
ytickformat(ax5, '%.1f');

f = fit(avgOverlap', avgRM', 'poly1');
plot(ax5.XLim, f(ax5.XLim), 'k')

% Save the remaining figures
fh1.Color = [1 1 1];
set(fh1, 'PaperPositionMode', 'auto')
set(fh1, 'InvertHardcopy', 'off')

fh2.Color = [1 1 1];
set(fh2, 'PaperPositionMode', 'auto')
set(fh2, 'InvertHardcopy', 'off')

fh3.Color = [1 1 1];
set(fh3, 'PaperPositionMode', 'auto')
set(fh3, 'InvertHardcopy', 'off')

fh5.Color = [1 1 1];
set(fh5, 'PaperPositionMode', 'auto')
set(fh5, 'InvertHardcopy', 'off')

print(fh1, '../../figures/Figure_2/fig2a.tif', '-dtiff')
print(fh2, '../../figures/Figure_2/fig2b.tif', '-dtiff')
print(fh3, '../../figures/Figure_2/fig2c.tif', '-dtiff')
print(fh5, '../../figures/Figure_S5/figS5a.tif', '-dtiff')
