close all;
clearvars;

% Load parameters
loadDefaultParameters
Beta = 25;

% Allocate array to store results
PP = nan(2, 250);

% Loop over random seed
parfor i = 1:size(PP, 2)
    
    % Set the Random seed
    rng(i)
    
    % Determine which case to simulate
    B = {1, 2, [1 2]};
    P = {[], 1};
    
    % Compute number of RM systems
    RM = unique([B{:} P{:}]);
    
    % Number of bacteria and phages
    nB = numel(B);
    nP = numel(P);
    
    % Starting condiiton
    x0 = Delta/(Eta*Beta)*[1/nB*ones(nB, 1); 10/(nB*nP)*ones(nB*nP, 1)]; % Everyone starts with a population of delta/(eta*beta)
    
    % Define omega and cost
    omega_0 = 10.^(lb + (ub-lb) * rand(size(RM)));
    cost    = 1-f*rand(size(RM));
    
    % Run dynamics
    [~, P_end] = DynamicalSystem('ExtendedModel', B, P, false, omega_0, cost, {Alpha, Beta, Eta, Delta, C, T}, x0, true);
    
    % Store the results
    PP(:, i) = P_end;
    
end

% Create figure
fh = figure('Resize', 'off');
fh.Position(3) = 400;
fh.Position(4) = 472.5;

ax1 = subplot(2,1,1); hold on; box on;
ax2 = subplot(2,1,2); hold on; box on;

ax1.Position = [0.22 0.60 0.73 0.37];
ax2.Position = [0.22 0.18 0.73 0.37];


ax1.FontSize = 20;
ax2.FontSize = 20;
ax1.LineWidth = 1;
ax2.LineWidth = 1;

histogram(ax1, PP(1, :) ./ sum(PP), 0:0.1:0.5, 'Normalization', 'Probability', 'FaceColor', [0.8 0.8 0.8]);
histogram(ax2, PP(2, :) ./ sum(PP), 0.5:0.1:1, 'Normalization', 'Probability', 'FaceColor', [0.8 0.8 0.8]);

ax1.XLim = [0 1];
ax2.XLim = [0 1];

ax1.XTick = [];

ax1.XTickLabel = {};

ax1.YLim = [0 1];
ax2.YLim = [0 1];

ax1.YTick = 0:0.5:1;
ax2.YTick = 0:0.5:1;

xtickformat(ax2, '%.1f');

ytickformat(ax1, '%.1f');
ytickformat(ax2, '%.1f');

xlabel(ax2, 'Fraction of population')

ylabel(ax1, 'pmf.')
ylabel(ax2, 'pmf.')

% Save the figures
if ~exist('../../figures/Figure_S10', 'dir')
    mkdir('../../figures/Figure_S10')
end

fh.Color = [1 1 1];
set(fh, 'PaperPositionMode', 'auto')
set(fh, 'InvertHardcopy', 'off')

print(fh, '../../figures/Figure_S10/figS10bc.tif', '-dtiff')