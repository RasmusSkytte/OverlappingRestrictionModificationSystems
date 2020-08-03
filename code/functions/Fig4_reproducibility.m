close all; clearvars;

% Start with default parameters
loadDefaultParameters

% Define the number of runs for for the figure
nFig4  = 6*6; % 7 different K, 6 repeats

fh1 = figure('Resize', 'off');
ax1 = subplot(2, 2, [1 3]); hold on; box on;
ax2 = subplot(2, 2, 4); hold on; box on;
ax3 = subplot(2, 2, 2); hold on; box on;

fh2 = figure('Resize', 'off'); clf; 
ax4 = axes; hold on; box on;


xlabel(ax1, '# Additions')
ylabel(ax1, '{\color{red}P/(10*C)}, B/C, {\color{blue} D/\beta}, {\color{orange}{\langle}#RM\rangle}')
ax1.XScale = 'log';
ax1.FontSize = 20;
ax1.LineWidth = 1;
ax1.YLim = [0 5];

pR = plot(ax1, nan, nan, 'LineWidth', 1.5, 'Color', [1 0.65 0], 'DisplayName', 'median(#RM)');
pP = plot(ax1, nan, nan, 'r', 'LineWidth', 1.5, 'DisplayName', 'P / C');
pB = plot(ax1, nan, nan, 'k', 'LineWidth', 1.5, 'DisplayName', 'B / C');
pD = plot(ax1, nan, nan, 'b', 'LineWidth', 1.5, 'DisplayName', 'D / \beta');
ax1.XLim  = [1 iterations];
ax1.XTick = logspace(0, 10, 6);

xlabel(ax2, '# RM');
ax2.XLim = [0 8];
ax2.XTick = [0 4 8];
ax2.YTick = [0 1];
ax2.YLim = [0 1];
ax2.FontSize = 20;
ax2.LineWidth = 1;

ax3.XTick = [];
ax3.YTick = [];
ax3.LineWidth = 1;


ax1.Position = [0.07 0.2 0.5 0.75];
ax2.Position = [0.65 0.2 0.3 0.2];
ax3.Position = [0.7 0.45 0.2 0.5];


cc = lines(7);

% Save the figures
drawnow;
fh1.Position = [10 80 1260 480];
fh2.Position = [488 342 560 420];

if ~exist('../../figures/Figure_4_seeds', 'dir')
    mkdir('../../figures/Figure_4_seeds')
end

fh1.Color = [1 1 1];
set(fh1, 'PaperPositionMode', 'auto')
set(fh1, 'InvertHardcopy', 'off')

% Allocate arrays
avgOverlap = [];
avgRM      = [];

% Loop over data
for k = 1:nFig4
    
    r = mod(k-1, 7) + 1; % Determine run
    s = ceil(k / 7) - 1; % Determine seed
    
    % Select run
    switch r
        case 2
            RM = 1:50;
        case 3
            RM = 1:100;
        case 4
            RM = 1:200;
        case 5
            RM = 1:400;
        case 6
            RM = 1:800;
    end
    
    if r == 1
        lname = sprintf('../data/Fig4/RM_inf_seed_%d.mat', s);
    else
        lname = sprintf('../data/Fig4/RM_%d_seed_%d.mat', RM(end), s);
    end
    
    % Load the data
    if exist(lname, 'file')
        load(lname)
    else 
        fprintf('%s is missing\n', lname)
        continue;
    end

    
    % Update figure
    pP.XData = 1:numel(phages);
    pP.YData = phages/(10*C);
    
    pB.XData = 1:numel(bacteria);
    pB.YData = bacteria/C;
    
    pD.XData = 1:numel(diversity);
    pD.YData = diversity/Beta;
    
    pR.XData = 1:numel(mRM);
    pR.YData = mRM;
    
    if numel(ax2.Children) > 0
        delete(ax2.Children(1));
    end
    if numel(ax3.Children) > 0
        delete(ax3.Children(1));
    end
    
    if r == 1
        % Plot the histogram
        h = histogram(ax2, nRM, 0.5:10.5, 'Normalization', 'Probability', 'FaceColor', 'b');
        
        % Reconstruct the implied network
        RM = 1:sum(nRM);
        
        B_sample = cell(numel(nRM), 1);
        for i = 1:numel(nRM)
            ind = randi(numel(RM), 1, nRM(i));
            B_sample{i} = ind;
            RM(ind) = [];
        end
        
        % Draw the graph at the end of the run
        [A_ij, A_RM] = computeNetwork(B_sample);
        G = graph(A_RM);
        plotNetwork(ax3, G)
        
        overlap = 0;
        
        title(ax3, sprintf('K = \\infty, seed = %d', s))
    else
        h = histogram(ax2, cellfun(@numel, B_samples{end}), 0.5:10.5, 'Normalization', 'Probability', 'FaceColor', 'b');
        
        % Plot the network
        [A_ij, A_RM] = computeNetwork(B_samples{end});
        G = graph(A_RM);
        plotNetwork(ax3, G)
        
        % Compute overlap
        overlap = (arrayfun(@(r) sum([B_samples{end}{:}]==r), unique([B_samples{end}{:}])) - 1)/(size(A_ij, 1) - 1);
        
        title(ax3, sprintf('K = %d, seed = %d, {\\langle}p{\\rangle} = %.3f', RM(end), s, mean(overlap)))
        
    end
    
    % Store average measures
    avgOverlap(end + 1) = mean(overlap);
    avgRM(end + 1)      = mean(sum(A_ij > 0, 2));
    
    
    saveas(fh1, sprintf('../../figures/Figure_4_seeds/r_%d_s_%d.png', r, s))
end


scatter(ax4, avgOverlap, avgRM, 'r', 'filled')
xlabel(ax4, '{\langle}p\rangle')
ylabel(ax4, '{\langle}#RM\rangle')
ax4.Position = [0.165 0.19 0.77 0.78];
ax4.FontSize = 20;
ax4.LineWidth = 1;
ytickformat(ax4, '%.1f');

% f = fit(avgOverlap', avgRM', 'poly1');
% plot(ax4, ax4.XLim, f(ax4.XLim), 'k')

% Save the remaining figures
fh2.Color = [1 1 1];
set(fh2, 'PaperPositionMode', 'auto')
set(fh2, 'InvertHardcopy', 'off')

print(fh2, '../../figures/Figure_SX/figSXb.tif', '-dtiff')
