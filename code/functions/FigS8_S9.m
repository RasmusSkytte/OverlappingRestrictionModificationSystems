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
if ~exist('../../figures/Figure_S8', 'dir')
    mkdir('../../figures/Figure_S8')
end
if ~exist('../../figures/Figure_S9', 'dir')
    mkdir('../../figures/Figure_S9')
end

% Prepare figures
fh1 = figure('Resize', 'off'); clf;
axNet = axes;

fh2 = figure('Resize', 'off'); clf; 
ax2 = axes; hold on; box on;

fh3 = figure('Resize', 'off'); clf; 
ax3 = axes; hold on; box on;

fh1.Color = [1 1 1];
set(fh1, 'PaperPositionMode', 'auto')
set(fh1, 'InvertHardcopy', 'off')

% Update figure size
drawnow;
fh1.Position(2) = 200;
fh1.Position(4) = fh1.Position(3);

fh2.Position = [488 342 560 420];
fh3.Position = [488 342 750 420];

% Ensure the data is sorted
[s, sortIndex] = sort(s);
S = S(sortIndex, :);
A_ij = A_ij(sortIndex, :);

% Select test bacteria
testBacteria = {'Bifidobacterium', 'Streptococcus', 'Paenibacillus', 'Escherichia'};

% Loop over subsets to generate figures
for d = 1:numel(unique(s))
    
    % Skip to the chosen networks
    if ~any(strcmp(testBacteria, uniqueLabels{3}{d}))
        continue
    end
    
    % Locate the start and size of the subset
    is = find(s == d, 1);
    ns = sum(s == d);
    
    % Extract the subset
    A_ij_s = A_ij(is:is+(ns-1), :);
    
    % Filter the subset for unique combinations
    A_ij_s = unique(A_ij_s, 'rows');
    
    % Label the adjacency matrix
    A_ij_s = A_ij_s .* repmat(1:size(A_ij_s, 2), size(A_ij_s, 1), 1);
    
    % Allocate array
    nBootstraps = 100;
    overlaps = nan(nBootstraps, 10);
    nRMs     = nan(nBootstraps, 10);
    
    % Draw 10 subsets of increasing size
    rng(60);
    for k = 1:10
        
        parfor r = 1:nBootstraps
            A_ij_ss = A_ij_s(randperm(size(A_ij_s, 1), round(size(A_ij_s, 1) * k / 10)), :);
            
            
            % Convert to list of RM format
            Bss = cell(size(A_ij_ss, 1), 1);
            for i = 1:size(A_ij_ss, 1)
                Bss{i} = A_ij_ss(i, A_ij_ss(i, :) > 0);
            end
            
            % Compute relative overlap between species
            A_ij_ss = A_ij_ss > 0;
            
            % Draw the graph at the end of the run
            [~, A_RM_ss] = computeNetwork(Bss);
            
            if isempty(A_RM_ss) % No RM systems present
                continue;
            end
            
            % Skip last (only 1 way to sample 100 percent)
            if k == 10 && r > 1
                continue;
            end
            
            % Compute overlap
            overlap_ss = (arrayfun(@(r) sum([Bss{:}]==r), unique([Bss{:}])) - 1)/(size(A_ij_ss, 1) - 1);
            overlaps(r, k) = mean(overlap_ss);
            
            % Store the average RM count
            nRMs(r, k) = mean(cellfun(@numel, Bss));
        end
    end
    
    % Create summary figures
    errorbar(ax2, (1:size(overlaps, 2)) * 10, nanmean(overlaps), nanstd(overlaps) ./ sqrt(sum(~isnan(overlaps))), 'o', 'MarkerFaceColor', 'auto', 'LineWidth', 2, 'DisplayName', uniqueLabels{g}{d})
    errorbar(ax3, (1:size(nRMs,     2)) * 10, nanmean(nRMs),     nanstd(nRMs)     ./ sqrt(sum(~isnan(nRMs))),     'o', 'MarkerFaceColor', 'auto', 'LineWidth', 2, 'DisplayName', uniqueLabels{g}{d})
    
    
    % Graph the most representitive network for each subsample frequency
    % for Escherichia
    if ~strcmp('Escherichia', uniqueLabels{3}{d})
        continue
    end
    
    % Draw 10 subsets of increasing size
    rng(60);
    for k = 1:10
        
        % Determine the most representative network
        [~, I] = min(overlaps(:, k) - nanmean(overlaps(:, k)));
        
        for r = 1:nBootstraps
            A_ij_ss = A_ij_s(randperm(size(A_ij_s, 1), round(size(A_ij_s, 1) * k / 10)), :);
            
            % Skip the remaning networks
            if r ~= I
                continue
            end
            
            % Convert to list of RM format
            Bss = cell(size(A_ij_ss, 1), 1);
            for i = 1:size(A_ij_ss, 1)
                Bss{i} = A_ij_ss(i, A_ij_ss(i, :) > 0);
            end
            
            % Compute relative overlap between species
            A_ij_ss = A_ij_ss > 0;
            
            % Draw the graph at the end of the run
            [~, A_RM_ss] = computeNetwork(Bss);
            
            if isempty(A_RM_ss) % No RM systems present
                continue;
            end
                        
            % Plot the RM network
            G = graph(A_RM_ss);
            l = arrayfun(@(d) num2str(d), find(sum(A_ij_ss, 1) > 0)', 'UniformOutput', false);
            l(G.degree < quantile(G.degree, 0.7)) = {' '};
            plotNetwork(axNet, G, l)

            % Adjust axes
            axNet.Position = [0 0 1 1];
            axNet.Visible = 'off';
            axNet.XLim = [-1.4 1.4];
            axNet.YLim = [-1.4 1.4]+0.2;

            % Label the figure
            text(axNet, -1.3,  1.50, sprintf('%s (%d percent)', labels{d}(2:end), 10*k), 'FontSize', 18, 'FontWeight', 'Bold');
            text(axNet, -1.3,  1.35, sprintf('# RM = %d,', size(A_RM_ss, 1)), 'FontSize', 14);
            text(axNet, -1.3,  1.25, sprintf('# samples = %d', size(A_ij_ss, 1)), 'FontSize', 14);
            text(axNet, -0.4,  1.35, sprintf('{\\langle}#RM{\\rangle} = %.2f,', mean(cellfun(@numel, Bss))), 'FontSize', 14);
            text(axNet, -0.4,  1.25, sprintf('{\\langle}p{\\rangle} = %.4f,',  overlaps(I, k)), 'FontSize', 14);
            
            % Save the network figure
            print(fh1, sprintf('../../figures/Figure_S9/figS9%s_%d_percent.tif', labels{d}, 10*k), '-dtiff')

            % Report overlap
            fprintf('<p_{%s}> = %.4f, n = %d\n', uniqueLabels{g}{d}, overlaps(I, k), size(A_ij_ss, 1))

        end
    end
end

xlabel(ax2, 'Percent sampled')
ylabel(ax2, '{\langle}p{\rangle}')
ax2.XLim = [0 100];
ax2.Position = [0.21 0.18 0.74 0.79];
ax2.FontSize = 20;
ax2.LineWidth = 1;
ytickformat(ax2, '%.3f');


l = legend(ax3, 'Location', 'NorthEastOutside');
xlabel(ax3, 'Percent sampled')
ylabel(ax3, '{\langle}#RM{\rangle}')
ax3.XLim = [0 100];
ax3.Position = [0.11 0.18 0.55 0.79];
ax3.FontSize = 20;
ax3.LineWidth = 1;
ytickformat(ax3, '%.1f');
l.Position(1) = 0.67;



fh2.Color = [1 1 1];
set(fh2, 'PaperPositionMode', 'auto')
set(fh2, 'InvertHardcopy', 'off')

fh3.Color = [1 1 1];
set(fh3, 'PaperPositionMode', 'auto')
set(fh3, 'InvertHardcopy', 'off')

print(fh2, '../../figures/Figure_S8/figS8a.tif', '-dtiff')
print(fh3, '../../figures/Figure_S8/figS8b.tif', '-dtiff')