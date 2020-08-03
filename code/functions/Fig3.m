close all;
clearvars;

% Set seed
rng(7);

% Load parameters
loadDefaultParameters
Beta = 25;

% Loop over cases to plot
labels = {'a', 'b', 'c'};
for i = 1:3
    % Determine which case to simulate
    switch i
        case 1
            B = {1, 2};
            P = {[]};
        case 2
            B = {1, 2, [1 2]};
            P = {[]};
        case 3
            B = {1, 2, [1 2]};
            P = {1};
    end

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
    [B_end, P_end, y, t] = DynamicalSystem('ExtendedModel', B, P, false, omega_0, cost, {Alpha, Beta, Eta, Delta, C, T}, x0, true);
    t(1) = eps;
    
    fh = figure('Resize', 'off'); 
    ax = gca; hold on; box on;
    ax.Position = [0.16 0.25 0.8 0.7];
    ax.FontSize = 20;
    ax.LineWidth = 1;
    ax.YTick = 10.^(0:2:8);
    cc = lines(3);
    cc(2, :) = [];
    for b = 1:nB
        plot(t, y(b, :), 'k', 'LineWidth', 4)
        if b <= 2
            plot(t, y(b, :), 'LineWidth', 3, 'Color', cc(b, :))
        else
            plot(t, y(b, :), 'LineWidth', 3, 'Color', cc(1, :))
            plot(t, y(b, :), '--', 'LineWidth', 3, 'Color', cc(2, :))
        end
    end
    set(gca, 'YScale', 'Log')
    set(gca, 'XScale', 'Log')
    xlim([1 t(end)])
    ax1.YLim(1) = 1;
    ylim([1 C])
    xlabel('Time (generations)')
    ylabel('Population Size')

    % Save the figures
    drawnow;
    fh.Position = [488 342 560 420];
    fh.Position(4) = 0.8 * fh.Position(4);
    
    if ~exist('../../figures/Figure_3', 'dir')
        mkdir('../../figures/Figure_3')
    end

    fh.Color = [1 1 1];
    set(fh, 'PaperPositionMode', 'auto')
    set(fh, 'InvertHardcopy', 'off')

    print(fh, sprintf('../../figures/Figure_3/fig3%s.tif', labels{i}), '-dtiff')

end
