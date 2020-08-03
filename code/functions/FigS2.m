close all;
clearvars;

% Load parameters
loadDefaultParameters
Beta = 25;

% Number of repititions
N = 17;

% Define test case
B = {1, 2, [1 2]};
P = {[]};

% Compute number of RM systems
RM = unique([B{:} P{:}]);

% Number of bacteria and phages
nB = numel(B);
nP = numel(P);

% Set up test space
[omega, gamma]   = meshgrid(10.^(0:-0.125:-2), linspace(0, 1, N));

% Prepare figure
fh = figure('Resize', 'off'); clf;
fh.Position = [0.5 0.2 2.25 1] .* fh.Position;

ax1 = subplot(1, 3, 1); hold on; box on;
ax2 = subplot(1, 3, 2); hold on; box on;
ax3 = subplot(1, 3, 3); hold on; box on;

ax1.Position = [0.08 0.16 0.25 0.83];
ax2.Position = [0.41 0.16 0.25 0.83];
ax3.Position = [0.74 0.16 0.25 0.83];

ax1.FontSize = 20;
ax2.FontSize = 20;
ax3.FontSize = 20;

ax1.LineWidth = 1;
ax2.LineWidth = 1;
ax3.LineWidth = 1;

ax1.TickLength = 1.5*ax1.TickLength;
ax2.TickLength = 1.5*ax2.TickLength;
ax3.TickLength = 1.5*ax3.TickLength;

ax1.YAxis.MinorTick = 'on';
ax1.YAxis.MinorTickValues = unique([log10(linspace(10^-2, 10^-1, 11)) log10(linspace(10^-1, 1, 11))]);

ax2.YAxis.MinorTick = 'on';
ax2.YAxis.MinorTickValues = unique([log10(linspace(10^-2, 10^-1, 11)) log10(linspace(10^-1, 1, 11))]);

ax3.YAxis.MinorTick = 'on';
ax3.YAxis.MinorTickValues = unique([log10(linspace(10^-2, 10^-1, 11)) log10(linspace(10^-1, 1, 11))]);

ax1.XLim = [-0.035 1.035];
ax1.YLim = [-2.0625 0.0625];

ax2.XLim = [-0.035 1.035];
ax2.YLim = [-2.0625 0.0625];

ax3.XLim = [-0.035 1.035];
ax3.YLim = [-2.0625 0.0625];

xlabel(ax1, '\gamma')
ylabel(ax1, '\omega')

xlabel(ax2, '\gamma')
xlabel(ax3, '\gamma')

ax1.XLabel.Position(2) = -2.25;
ax2.XLabel.Position(2) = -2.25;
ax3.XLabel.Position(2) = -2.25;

% Prepare data arrays
ps = zeros(size(omega));
Biomass   = zeros(size(omega));
B_double  = zeros(size(omega));

% Loop over omega, gamma
indicies = randperm(numel(omega));
for i = indicies(1:end)

    % Load data
    lname = sprintf('../data/FigS2/omega_1e%.3f_gamma_%.4f.mat', log10(omega(i)), gamma(i));
    if exist(lname, 'file')
        load(sprintf('../data/FigS2/omega_1e%.3f_gamma_%.4f.mat', log10(omega(i)), gamma(i)));
    else
        continue
    end

    % Set booleans
    b_AB = false;
    b_A_B = false;
    b_A_B_AB = false;
    
    % Set up starting coordinates
    deb = Delta/(Eta*(Beta-1));        
    deb2 = 1.9*Delta / (Eta*(Beta*(1+omega(i))-2)); % Use 1.9 instead of 2 to break degeneracy when omega  = 1
    b = [deb deb2 1e1 1e3 1e5 1e7];

    [pi, b1, b2] = meshgrid(10*b, b, b);

    % Loop over starting conditions
    for j = 1:numel(pi)

        % Retrieve B_end and P_end
        B_end = BB(and(and(coordinates(:, 1)==b1(j), coordinates(:, 2)==b2(j)), coordinates(:, 3)==pi(j)), :);
        P_end = PP(and(and(coordinates(:, 1)==b1(j), coordinates(:, 2)==b2(j)), coordinates(:, 3)==pi(j)));
        Biomass(i) = sum(B_end);
        B_double(i) = B_end(end);

        % Determine posible soluitions
        thres = 1;
        if B_end(1) > thres && B_end(end) < thres
            if ~b_A_B
                b_A_B = true;
                ps(i) = ps(i) + 2^0;
            end
        end

        if B_end(1) > thres && B_end(end) > thres
            if ~b_A_B_AB
                b_A_B_AB = true;
                ps(i) = ps(i) + 2^1;
            end
        end

        if B_end(1) < thres && B_end(end) > thres
            if ~b_AB
                b_AB = true;
                ps(i) = ps(i) + 2^2;
            end
        end
    end
end

% Plot the phasespace
% Test if points have reached SS
if any(ps(:) == 3)
    error('SS not reached')
end
if any(ps(gamma==Alpha) ~= 0)
    error('SS not reached')
end

% Relabel the points
ps(ps == 6) = 3;
ps(ps == 7) = 4;

% Plot solution space
imagesc(ax1, gamma(:, 1), log10(omega(1, :)), ps');
colormap(ax1, [0.7 0.7 0.7; lines(4)])
ax1.CLim = [0 4];

ax0 = axes;
ax0.Visible = 'off';
ax0.Position = [0 0 1 1];

text(ax0, 0.11,  0.575, 'No Solution', 'Color', 'w', 'FontWeight', 'Bold', 'FontSize', 18, 'HorizontalAlignment', 'Center', 'Rotation', 90);
text(ax0, 0.175, 0.87,  'A/B',         'Color', 'w', 'FontWeight', 'Bold', 'FontSize', 18, 'HorizontalAlignment', 'Center');
text(ax0, 0.255, 0.40,  'A/B/AB',      'Color', 'w', 'FontWeight', 'Bold', 'FontSize', 18, 'HorizontalAlignment', 'Center');
text(ax0, 0.32,  0.575, 'AB',          'Color', 'w', 'FontWeight', 'Bold', 'FontSize', 18, 'HorizontalAlignment', 'Center', 'Rotation', 90);

% Plot change in biomass
B_ref = 2*Delta ./ (Eta*(Beta*(1+omega)-2));
DeltaB = (Biomass-B_ref)./B_ref;
imagesc(ax2, gamma(5:end, 1), log10(omega(1, :)), DeltaB(6:end, :)');
ax2.Color = [0.7 0.7 0.7];
cc = spring(128);
colormap(ax2, cc(33:96, :))

% Plot relative biomass of AB solution
RelB = B_double./Biomass;
imagesc(ax3, gamma(5:end, 1), log10(omega(1, :)), RelB(6:end, :)');
ax3.Color = [0.7 0.7 0.7];
colormap(ax3, cc(65:128, :))

% Fix the ticks
ax1.YTick = -2:0;
ax1.YTickLabel = cellfun(@(s)sprintf('10^{%s}', s), ax1.YTickLabel, 'UniformOutput', false);
set(ax1, 'Layer', 'top')

ax2.YTick = -2:0;
ax2.YTickLabel = cellfun(@(s)sprintf('10^{%s}', s), ax2.YTickLabel, 'UniformOutput', false);
set(ax2, 'Layer', 'top')

ax3.YTick = -2:0;
ax3.YTickLabel = cellfun(@(s)sprintf('10^{%s}', s), ax3.YTickLabel, 'UniformOutput', false);
set(ax3, 'Layer', 'top')



% Guide lines
k = @(o) 1 - 2*Delta/(C*Eta*(Beta*(1+o)-2));
D = @(o) o^2*k(o)^2 + k(o)^2*Alpha*(1+o)*(1-o);
f = @(o) (o*k(o)+sqrt(D(o)))/((1+o)*k(o));
ww = logspace(-2, 0);
gg = arrayfun(f, ww);

plot(ax1, [1 1]*Alpha, ax1.YLim, 'w-.', 'LineWidth', 2, 'DisplayName', '$\gamma=\alpha$');
plot(ax2, [1 1]*Alpha, ax2.YLim, 'w-.', 'LineWidth', 2, 'DisplayName', '$\gamma=\alpha$');
plot(ax3, [1 1]*Alpha, ax3.YLim, 'w-.', 'LineWidth', 2, 'DisplayName', '$\gamma=\alpha$');

plot(ax1, gg, log10(ww), 'w-', 'LineWidth', 2, 'DisplayName', '$\gamma = f(\omega)$');
plot(ax2, gg, log10(ww), 'w-', 'LineWidth', 2, 'DisplayName', '$\gamma = f(\omega)$');
plot(ax3, gg, log10(ww), 'w-', 'LineWidth', 2, 'DisplayName', '$\gamma = f(\omega)$');

c1 = colorbar(ax2, 'Location', 'West');
ax2.CLim(1) = -0.5;
ax2.CLim(2) = 0.5;
c1.TickLabels{2} = '0.0';
c1.Position(1) = c1.Position(1)-0.0034;
c1.Position(3) = 0.01;
c1.Color = 'w';

c2 = colorbar(ax3, 'Location', 'West');
ax3.CLim(2) = 1;
c2.Ticks = 0:0.5:1;
c2.TickLabels{1} = '0.0';
c2.TickLabels{end} = '1.0';
c2.Position(1) = c2.Position(1)-0.0034;
c2.Position(3) = 0.01;
c2.Color = 'w';

% Save the figures
if ~exist('../../figures/Figure_S2', 'dir')
    mkdir('../../figures/Figure_S2')
end

fh.Color = [1 1 1];
set(fh, 'PaperPositionMode', 'auto')
set(fh, 'InvertHardcopy', 'off')

print(fh, '../../figures/Figure_S2/figS2.tif', '-dtiff')