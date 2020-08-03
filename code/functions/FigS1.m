close all;
clearvars;

% Set seed
rng(60);

% Number of repititions
M = 30;

% Test #
test = 2;

% Load parameters
loadDefaultParameters
Beta = 25;
T = 1e4;

% Compute SS value
deb = Delta/(Eta*(Beta-1));

% Define simple case
switch test
    case 1
        B = {1, 2};
        P = {[]};
    case 2
        B = {1, 2, [1 2]};
        P = {[]};
end

% Define test omegas
omegas = 10.^[-0.25 -0.5 -1 -2];

% Compute number of RM systems
RM = unique([B{:} P{:}]);

% Number of bacteria and phages
nB = numel(B);
nP = numel(P);

% Prepare figure
fh = figure(1); clf;
fh.Position = [0.5 0.2 2.25 1.5] .* fh.Position;
ax1 = subplot(2, 2, [1 3]); hold on; box on;

xlabel(ax1, '\gamma')
ylabel(ax1, 'B')

ax2 = subplot(2, 2, 2); hold on; box on;
xlabel(ax2, '\gamma')
ylabel(ax2, 'b_A, b_B')

ax3 = subplot(2, 2, 4); hold on; box on;
xlabel(ax3, '\gamma')
ylabel(ax3, 'b_{AB}')

ax1.Position = [0.08 0.13  0.50 0.815];
ax2.Position = [0.66 0.625 0.33 0.3175];
ax3.Position = [0.66 0.13  0.33 0.3175];

ax1.FontSize = 20;
ax2.FontSize = 20;
ax3.FontSize = 20;

ax1.LineWidth = 1;
ax2.LineWidth = 1;
ax3.LineWidth = 1;

ax1.XLim = [0 1];
ax2.XLim = [0 1];
ax3.XLim = [0 1];

ax1.YLim = [0 round((3*deb)/10^floor(log10(3*deb)), 2) * 10^floor(log10(3*deb))];
ax2.YLim = [0  ceil((  deb)/10^floor(log10(  deb)))    * 10^floor(log10(  deb))];
ax3.YLim = [0  ceil((  deb)/10^floor(log10(  deb)))    * 10^floor(log10(  deb))];

% text(ax1, 0, 2.6e6, '{\times}10^6', 'FontSize', 20)

% Symbols
symbols = {'v', 's', 'o', '^'};

% Guide lines
plot(ax1, [1 1]*Alpha, ax1.YLim, 'k-.', 'LineWidth', 1.5);
plot(ax2, [1 1]*Alpha, ax2.YLim, 'k-.', 'LineWidth', 1.5, 'DisplayName', '$\gamma = \alpha$');
plot(ax3, [1 1]*Alpha, ax3.YLim, 'k-.', 'LineWidth', 1.5);

f = @(g, o)(g^2*(1+o)-2*o*g)*(1-2*Delta/(Beta*(1+o)-2))-Alpha*(1-o);
gg = arrayfun(@(o)fzero(@(g)f(g, o), 1), omegas);

plot(ax1, reshape(repmat(gg, 2, 1), 8, 1), ax1.YLim(2)*repmat([-0.1 1.1 1.1 -0.1]', 2, 1), 'k--', 'LineWidth', 1.5);
plot(ax2, reshape(repmat(gg, 2, 1), 8, 1), ax2.YLim(2)*repmat([-0.1 1.1 1.1 -0.1]', 2, 1), 'k--', 'LineWidth', 1.5, 'DisplayName', '$\gamma = f(\omega)$');
plot(ax3, reshape(repmat(gg, 2, 1), 8, 1), ax3.YLim(2)*repmat([-0.1 1.1 1.1 -0.1]', 2, 1), 'k--', 'LineWidth', 1.5);


plot(ax1, ax1.XLim, [1 1]*deb,    'Color', 'k', 'LineWidth', 1.5);
plot(ax2, ax2.XLim, 10*[1 1]*deb, 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', '$\frac{\delta}{\eta(\beta-1)}$');
plot(ax3, ax3.XLim, [1 1]*deb,    'Color', 'k', 'LineWidth', 1.5);

plot(ax1, repmat([-0.1 1.1 1.1 -0.1]', 2, 1), 2*Delta ./ (Eta*(Beta*(1+reshape([omegas; omegas], 2*numel(omegas), 1))-2)), 'k:', 'LineWidth', 1.5);
plot(ax2, repmat([-0.1 1.1 1.1 -0.1]', 2, 1), Delta ./ (Eta*(Beta*(1+reshape([omegas; omegas], 2*numel(omegas), 1))-2)), 'k:', 'LineWidth', 1.5, 'DisplayName', '$\frac{2\delta}{\eta(\beta(1+\omega)-2)}$');

cc = lines(numel(omegas));

% Loop over omegas
for o = omegas

    s1 = scatter(ax1, nan, nan, symbols{1}, 'MarkerEdgeColor', cc(1, :), 'MarkerFaceColor', cc(1, :), 'LineWidth', 2, 'DisplayName', ['\omega = 10^{' num2str(log10(o)) '}']);
    s2 = scatter(ax2, nan, nan, symbols{1}, 'MarkerEdgeColor', cc(1, :), 'MarkerFaceColor', cc(1, :), 'LineWidth', 2);
    s3 = scatter(ax3, nan, nan, symbols{1}, 'MarkerEdgeColor', cc(1, :), 'MarkerFaceColor', cc(1, :), 'LineWidth', 2);
    symbols(1) = [];
    cc(1, :) = [];

    % Loop over gammas
    for r = linspace(0, 1, M)

        % Loop over starting conditions
        for s = 1

            % Starting condiiton
            switch s
                case 1
                    % delta eta beta limit
                    x0 = deb*[1/nB*ones(nB, 1); 10/(nB*nP)*ones(nB*nP, 1)];
                case 2
                    % single member limit
                    x0 = [ones(nB, 1); ones(nB*nP, 1)];
                case 3
                    % max biomass limit
                    x0 = [1/nB*ones(nB, 1)*C; ones(nB*nP, 1)];
                otherwise
                    error('Starting condiiton not known!')
            end


            % Define omega and cost
            omega_0 = ones(size(RM)) * o;
            cost    = ones(size(RM)) * r;

            % Run dynamics
            B_end = DynamicalSystem('ExtendedModel', B, P, false, omega_0, cost, {Alpha, Beta, Eta, Delta, C, T}, x0, true);
            
            % Plot the outcomes
            s1.XData = [s1.XData cost(1)];
            s1.YData = [s1.YData sum(B_end)];

            s2.XData = [s2.XData cost];
            s2.YData = [s2.YData B_end(1:2)'];

            s3.XData = [s3.XData cost(1)];
            s3.YData = [s3.YData B_end(3)];

            drawnow;
        end
    end
end
l1 = legend(ax1, ax1.Children(1:4),       'Location', 'NorthEastOutside', 'FontSize', 14);
l2 = legend(ax2, ax2.Children(end-3:end), 'Location', 'NorthEastOutside', 'FontSize', 14);

l1.Position = [0.44 0.725 0.13 0.19];
l2.Position = [0.42 0.148 0.15 0.23];
l2.Interpreter = 'Latex';
l2.FontSize = 20;

% Save the figures
if ~exist('../../figures/Figure_S1', 'dir')
    mkdir('../../figures/Figure_S1')
end

fh.Color = [1 1 1];
set(fh, 'PaperPositionMode', 'auto')
set(fh, 'InvertHardcopy', 'off')

print(fh, '../../figures/Figure_S1/figS1.tif', '-dtiff')
