function [B_end, P_end, D_end, bacteria, phages, diversity, mRM, gamma, omega, B_samples, nRM, age] = simulateOriginalModel(Alpha, Beta, Eta, Delta, C, T, S, f, lb, ub, iterations, fname, varargin)

% Initialize empty lists
gamma = [];
omega = [];
nRM   = [];
age   = [];

I = [];
y = [];

bacteria  = nan(1, iterations);
phages    = nan(1, iterations);
diversity = nan(1, iterations);
mRM       = nan(1, iterations);

if numel(varargin) == 2
    pertubation = varargin{1};
    sampleTimes = varargin{2};
else
    pertubation = nan;
    sampleTimes = [];
end

B_samples = {};

% Main loop
for i = 1:iterations

    % Generate new species
    if i <= iterations

        % Generate new gamma, omega value
        [g, o, n] = generateSpecies(nRM, lb, ub, f, S, 1);

        gamma = [gamma; g'];
        omega = [omega; o'];
        nRM   = [nRM;   n'];

        % Store the age
        age = [age; zeros(numel(gamma)-numel(age), 1)];

    end

    % Count the number of species
    nB = numel(gamma);

    % Set starting conditions
    if i == 1
        x0 = ones(2*nB, 1); % Everyone starts with a population of 1
    else
        b = ones(nB-sum(I), 1); % New bacteria and phages start with a population of 1
        x0 = [y(1:sum(I), end); b; abs(y(sum(I)+(1:sum(I)), end)); b];
    end

    % Run dynamics
    [~, ~, y] = DynamicalSystem('OriginalModel', {[]}, {[]}, false, omega, gamma, {Alpha, Beta, Eta, Delta, C, T}, x0);

    % Determine the end populations
    B_end = y(1:nB, end);
    P_end = y((nB+1):end, end);

    % Apply threshold of survival
    I = B_end > 1;

    % Ignore species below threshold
    B_end(~I) = [];

    % Store data
    diversity(i) = sum(I);
    bacteria(i)  = sum(B_end);
    phages(i)    = sum(P_end);
    D_end        = diversity(i);
    mRM(i)       = mean(nRM);

    % Remove the dead species
    y([~I; ~I], :) = [];
    gamma(~I) = [];
    omega(~I) = [];
    nRM(~I) = [];
    age(~I) = [];

    % Increase the age
    age = age + 1;

    % Sample B if required
    if any(i == sampleTimes)
        B_samples{end+1} = nRM;
    end

    if mod(i, 100) == 0
        save(fname, 'bacteria', 'phages', 'diversity', 'mRM', 'nRM', 'gamma', 'omega', 'age', 'B_samples', 'C', 'Alpha', 'Beta', 'Delta', 'T', 'lb', 'ub', 'S', 'f', 'iterations', 'sampleTimes')
    end
end
