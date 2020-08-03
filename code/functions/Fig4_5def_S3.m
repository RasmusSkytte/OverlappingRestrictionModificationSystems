close all;
clearvars;

% Load parameters
loadDefaultParameters

% Prepare figure
close all;
fh1 = figure('Resize', 'off'); clf;
ax1 = subplot(2, 3, 1); hold on; box on;
ax2 = subplot(2, 3, 2); hold on; box on;
ax3 = subplot(2, 3, 3); hold on; box on;
ax4 = subplot(2, 3, 4); hold on; box on;
ax5 = subplot(2, 3, 5); hold on; box on;
ax6 = subplot(2, 3, 6); hold on; box on;

ax1.Position = [0.06  0.48 0.27 0.5];
ax2.Position = [0.385 0.48 0.27 0.5];
ax3.Position = [0.71  0.48 0.27 0.5];
ax4.Position = [0.06  0.13 0.27 0.2];
ax5.Position = [0.385 0.13 0.27 0.2];
ax6.Position = [0.71  0.13 0.27 0.2];

xlabel(ax1, '# Additions');
xlabel(ax2, '# Additions');
xlabel(ax3, '# Additions');
xlabel(ax4, '# RM');
xlabel(ax5, '# RM');
xlabel(ax6, '# RM');

ax0 = axes;
ax0.Visible = 'off';
ax0.Position = [0 0 1 1];

text(ax0, 0.015, 0.70, '{\color{red}P/(10*C)}, B/C, {\color{blue} D/\beta}', 'FontSize', 18, 'Rotation', 90, 'HorizontalAlignment', 'Center')
text(ax0, 0.015, 0.22, 'pmf', 'FontSize', 18, 'Rotation', 90, 'HorizontalAlignment', 'Center')

ax1.XScale = 'log';
ax2.XScale = 'log';
ax3.XScale = 'log';

ax1.FontSize = 20;
ax2.FontSize = 20;
ax3.FontSize = 20;
ax4.FontSize = 20;
ax5.FontSize = 20;
ax6.FontSize = 20;

ax1.LineWidth = 1;
ax2.LineWidth = 1;
ax3.LineWidth = 1;
ax4.LineWidth = 1;
ax5.LineWidth = 1;
ax6.LineWidth = 1;

ax4.XLim = [0 8];
ax5.XLim = [0 8];
ax6.XLim = [0 8];

ax4.XTick = [0 4 8];
ax5.XTick = [0 4 8];
ax6.XTick = [0 4 8];

ax1.YLim = [0 3];
ax2.YLim = [0 3];
ax3.YLim = [0 3];
ax4.YLim = [0 1];
ax5.YLim = [0 1];
ax6.YLim = [0 1];

ax1.YTick = 0:5;
ax2.YTick = 0:5;
ax3.YTick = 0:5;
ax4.YTick = 0:1;
ax5.YTick = 0:1;
ax6.YTick = 0:1;

fh2 = figure('Resize', 'off');

ax8 = axes;

ax8.Position = [0 0 1 1];

fh3 = figure('Resize', 'off'); 

ax9 = axes;

ax9.Position = [0 0 1 1];

fh4 = figure('Resize', 'off');

ax10 = axes;

ax10.Position = [0 0 1 1];


% Plot the histogram of overlap
fh5 = figure('Resize', 'off');
ax11 = axes; hold on; box on;

xlabel(ax11, '# Intersection / # Union')
ylabel(ax11, 'pmf')

ax11.Position = [0.155 0.18 0.84 0.79];
ax11.FontSize = 20;
ax11.LineWidth = 1;

% Plot the histogram of min overlap
fh6 = figure('Resize', 'off');
ax12 = axes; hold on; box on;

xlabel(ax12, 'max(# Intersection / # Union)')
ylabel(ax12, 'pmf')

ax12.Position = [0.155 0.18 0.84 0.79];
ax12.FontSize = 20;
ax12.LineWidth = 1;


for j = 1:3

    % Set settings for the run
    switch j
        case 1
            axPop  = ax1;
            axHist = ax4;
            axNet  = ax8;
        case 2
            RM = 1:800;
            s = 4;
            axPop  = ax2;
            axHist = ax5;
            axNet  = ax9;
        case 3
            RM = 1:50;
            s = 5;
            axPop  = ax3;
            axHist = ax6;
            axNet  = ax10;
    end

    % Prepare figure
    pP = plot(axPop, nan, nan, 'r', 'LineWidth', 1.5, 'DisplayName', 'P / C');
    pB = plot(axPop, nan, nan, 'k', 'LineWidth', 1.5, 'DisplayName', 'B / C');
    pD = plot(axPop, nan, nan, 'b', 'LineWidth', 1.5, 'DisplayName', 'D / \beta');

    axPop.XLim  = [1 iterations];
    axPop.XTick = logspace(0, 10, 6);

    if j == 1

        lname = '../data/Fig4/RM_inf_seed_0.mat';
        if exist(lname, 'file')
            load(lname)

            % Update figure
            pP.XData = 1:numel(phages);
            pP.YData = phages/(10*C);

            pB.XData = 1:numel(bacteria);
            pB.YData = bacteria/C;

            pD.XData = 1:numel(diversity);
            pD.YData = diversity/Beta;

        else 
            fprintf('%s is missing\n', lname)
            continue;
        end

        % Plot the histogram
        h = histogram(axHist, nRM, 0.5:10.5, 'Normalization', 'Probability');
        
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
        plotNetwork(axNet, G) 

    else

        lname = sprintf('../data/Fig4/RM_%d_seed_%d.mat', RM(end), s);
        if exist(lname, 'file')
            load(lname)

            % Update figure
            pP.XData = 1:numel(phages);
            pP.YData = phages/(10*C);

            pB.XData = 1:numel(bacteria);
            pB.YData = bacteria/C;

            pD.XData = 1:numel(diversity);
            pD.YData = diversity/Beta;

        else 
            fprintf('%s is missing\n', lname)
            continue;
        end

        % Plot the histogram
        histogram(axHist, cellfun(@numel, B_samples{end}), 0.5:10.5, 'Normalization', 'Probability');

        % Draw the graph at the end of the run
        [A_ij, A_RM] = computeNetwork(B_samples{end});
        G = graph(A_RM);
        plotNetwork(axNet, G)
        
        % Plot the "fig2-like" histograms
        % Compute overlap
        overlap = (arrayfun(@(r) sum([B_samples{end}{:}]==r), unique([B_samples{end}{:}])) - 1)/(size(A_ij, 1) - 1);

        % Report overlap
        fprintf('<p_RM_%d> = %.4f\n', RM(end), mean(overlap))

        % Compute relative overlap between species
        o_ij = nan(numel(B_samples{end}), numel(B_samples{end}));
        for k = 1:numel(B_samples{end})
            for l = 1:numel(B_samples{end})
                if k == l
                    o_ij(k, l) = nan;
                    continue;
                end
                o_ij(k, l) = sum(and(A_ij(k, :), A_ij(l, :))) / sum( or(A_ij(k, :), A_ij(l, :)));
            end
        end

       
        % Plot the histogram of overlap
        histogram(ax11, o_ij(:), 0:0.2:1, 'LineWidth', 1, 'FaceAlpha', 0.4, 'Normalization', 'Probability');

        % Plot the histogram of min overlap
        histogram(ax12, max(o_ij), 0:0.2:1, 'LineWidth', 1, 'FaceAlpha', 0.4, 'Normalization', 'Probability');
        
    end
    
    
    % Adjust axes
    axNet.Position = [0 0 1 1];
    axNet.Visible = 'off';
    axNet.XLim = [-1.4 1.4];
    axNet.YLim = [-1.4 1.4]+0.2;

    % Label the figure
    if j == 1
        text(axNet, -1.3,  1.50, 'K \rightarrow \infty', 'FontSize', 18, 'FontWeight', 'Bold');
        text(axNet, -0.4,  1.25, '{\langle}p{\rangle} = 0', 'FontSize', 14);
    else
        text(axNet, -1.3,  1.50, sprintf('K = %d', RM(end)), 'FontSize', 18, 'FontWeight', 'Bold');
        text(axNet, -0.4,  1.25, sprintf('{\\langle}p{\\rangle} = %.4f,', mean(overlap)), 'FontSize', 14);
    end
    text(axNet, -1.3, 1.35, sprintf('# RM = %d,', size(A_RM, 1)), 'FontSize', 14);
    text(axNet, -0.4, 1.35, sprintf('{\\langle}#RM{\\rangle} = %.2f,', mean(sum(A_ij > 0, 2))), 'FontSize', 14);
    text(axNet, -1.3, 1.25, sprintf('# samples = %d', size(A_ij, 1)), 'FontSize', 14);

end


% Save the figures
drawnow;

fh1.Position = [10 80 1260 613];

fh2.Position(2) = 200;
fh2.Position(4) = fh2.Position(3);

fh3.Position(2) = 200;
fh3.Position(4) = fh2.Position(3);

fh4.Position(2) = 200;
fh4.Position(4) = fh2.Position(3);

fh5.Position = [488 342 560 420];
fh6.Position = [488 342 560 420];


if ~exist('../../figures/Figure_4', 'dir')
    mkdir('../../figures/Figure_4')
end

fh1.Color = [1 1 1];
set(fh1, 'PaperPositionMode', 'auto')
set(fh1, 'InvertHardcopy', 'off')

print(fh1, '../../figures/Figure_4/fig4.tif', '-dtiff')


if ~exist('../../figures/Figure_5', 'dir')
    mkdir('../../figures/Figure_5')
end

fh2.Color = [1 1 1];
set(fh2, 'PaperPositionMode', 'auto')
set(fh2, 'InvertHardcopy', 'off')

fh3.Color = [1 1 1];
set(fh3, 'PaperPositionMode', 'auto')
set(fh3, 'InvertHardcopy', 'off')

fh4.Color = [1 1 1];
set(fh4, 'PaperPositionMode', 'auto')
set(fh4, 'InvertHardcopy', 'off')

print(fh2, '../../figures/Figure_5/fig5d.tif', '-dtiff')
print(fh3, '../../figures/Figure_5/fig5e.tif', '-dtiff')
print(fh4, '../../figures/Figure_5/fig5f.tif', '-dtiff')


if ~exist('../../figures/Figure_S3', 'dir')
    mkdir('../../figures/Figure_S3')
end

fh5.Color = [1 1 1];
set(fh5, 'PaperPositionMode', 'auto')
set(fh5, 'InvertHardcopy', 'off')

fh6.Color = [1 1 1];
set(fh6, 'PaperPositionMode', 'auto')
set(fh6, 'InvertHardcopy', 'off')

print(fh5, '../../figures/Figure_S3/figS3a.tif', '-dtiff')
print(fh6, '../../figures/Figure_S3/figS3b.tif', '-dtiff')