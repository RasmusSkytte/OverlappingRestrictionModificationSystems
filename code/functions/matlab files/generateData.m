function generateData(k)

% Start with default parameters
loadDefaultParameters

% Make ../data folders
if ~exist(('../data/Fig6', 'dir')
    mkdir('../data/Fig6')
end

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
    fname = sprintf('../data/Fig6/RM_inf_seed_%d.mat', s);
    if ~exist(fname, 'file')
        [~ , ~, ~, bacteria, phages, diversity, mRM, gamma, omega, B_samples, nRM, age] = simulateOriginalModel(Alpha, Beta, Eta, Delta, C, T, S, f, lb, ub, iterations, fname, nan, sampleTimes);
        save(fname, 'bacteria', 'phages', 'diversity', 'mRM', 'nRM', 'gamma', 'omega', 'B_samples', 'age', 'C', 'Alpha', 'Beta', 'Delta', 'T', 'lb', 'ub', 'S', 'f', 'iterations', 'sampleTimes')
    end
else
    fname = sprintf('../data/Fig6/RM_%d_seed_%d.mat', RM(end), s);
    if ~exist(fname, 'file')
        [B, ~, ~, ~, ~, bacteria, phages, diversity, mRM, cost, omega_0, B_samples, overlap, age] = simulateExtendedModel(Alpha, Beta, Eta, Delta, C, T, S, RM, f, lb, ub, iterations, fname, nan, sampleTimes);
        save(fname, 'B', 'bacteria', 'phages', 'diversity', 'overlap', 'mRM', 'cost', 'omega_0', 'B_samples', 'age', 'C', 'Alpha', 'Beta', 'Delta', 'T', 'lb', 'ub', 'S', 'f', 'iterations', 'sampleTimes')
    end
end


end
