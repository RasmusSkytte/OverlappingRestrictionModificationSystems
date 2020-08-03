function generateData(k)

% Start with default parameters
loadDefaultParameters

% Make ../data folders
if ~exist(sprintf('../data/Fig4_Beta_%d', Beta), 'dir')
    mkdir(sprintf('../data/Fig4_Beta_%d', Beta))
end

if ~exist('../data/FigS2', 'dir')
    mkdir('../data/FigS2')
end

% Define the number of runs for each generator
nFig4  = 6*6; % 6 different K, 6 repeats
nFigS2 = 289; % 17 x 17 phasespace

% Check which scenario to generate ../data for
if k <= nFig4
        
    r = mod(k-1, 7) + 1; % Determine run
    s = ceil(k / 7) - 1; % Determine seed
    
    rng(s);     % Set seed

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

    sampleTimes = 10.^(0:log10(iterations));

    if r == 1
        sname = sprintf('../data/Fig4_Beta_%d/RM_inf_seed_%d.mat', Beta, s);
        if ~exist(sname, 'file')
            [~ , ~, ~, bacteria, phages, diversity, mRM, gamma, omega, B_samples, nRM, age] = simulateOriginalModel(Alpha, Beta, Eta, Delta, C, T, S, f, lb, ub, iterations, sname, nan, sampleTimes);
            save(sname, 'bacteria', 'phages', 'diversity', 'mRM', 'nRM', 'gamma', 'omega', 'B_samples', 'age', 'C', 'Alpha', 'Beta', 'Delta', 'T', 'lb', 'ub', 'S', 'f', 'iterations', 'sampleTimes')
        end
    else
        sname = sprintf('../data/Fig4_Beta_%d/RM_%d_seed_%d.mat', Beta, RM(end), s);
        if ~exist(sname, 'file')
            [B, ~, ~, ~, ~, bacteria, phages, diversity, mRM, cost, omega_0, B_samples, overlap, age] = simulateExtendedModel(Alpha, Beta, Eta, Delta, C, T, S, RM, f, lb, ub, iterations, sname, nan, sampleTimes);
            save(sname, 'B', 'bacteria', 'phages', 'diversity', 'overlap', 'mRM', 'cost', 'omega_0', 'B_samples', 'age', 'C', 'Alpha', 'Beta', 'Delta', 'T', 'lb', 'ub', 'S', 'f', 'iterations', 'sampleTimes')
        end
    end

elseif k > nFig4 && k <= (nFig4 + nFigS2)

    i = k - nFig4;

    N = 17;
    [omega, gamma]   = meshgrid(10.^(0:-0.125:-2), linspace(0, 1, N));

    % Define gamma and omega
    omega = omega(i);
    gamma = gamma(i);
    
    % Overwrite the beta value
    Beta = 25;

    fname = sprintf('../data/FigS2/omega_1e%.3f_gamma_%.4f.mat', log10(omega), gamma);
    if ~exist(fname, 'file')

        % Define test case
        B = {1, 2, [1 2]};
        P = {[]};

        % Compute number of RM systems
        RM = unique([B{:} P{:}]);

        % Number of bacteria and phages
        nB = numel(B);

        % Set up test space
        deb = Delta/(Eta*(Beta-1));
        deb2 = 1.9*Delta / (Eta*(Beta*(1+omega)-2));    % Use 1.9 instead of 2 to break degeneracy when omega  = 1
        b = [deb deb2 1e1 1e3 1e5 1e7];

        [pi, b1, b2] = meshgrid(10*b, b, b);

        % Allocate arrays
        coordinates = nan(numel(pi), 3);
        BB = nan(numel(pi), 3);
        PP = nan(size(pi));

        % Loop over starting conditions
        for j = 1:numel(pi)

            % Set starting conditions
            x0 = [b1(j); b1(j); b2(j); ones(nB, 1)*pi(j)];

            % Define omega and cost
            omega_0 = ones(size(RM)) * omega;
            cost    = ones(size(RM)) * gamma;

            % Run dynamics
            [B_end, P_end] = DynamicalSystem('ExtendedModel', B, P, false, omega_0, cost, {Alpha, Beta, Eta, Delta, C, 1e5}, x0, true);

            coordinates(j, :) = [b1(j) b2(j) pi(j)];
            BB(j, :) = B_end';
            PP(j) = P_end;

            % Save results
            save(fname, 'coordinates', 'BB', 'PP');

        end
    end
end
    
end
