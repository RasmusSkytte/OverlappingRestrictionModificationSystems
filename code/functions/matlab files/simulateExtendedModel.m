function [B, P, B_end, P_end, D_end, bacteria, phages, diversity, mRM, cost, omega_0, B_samples, overlap, age] = simulateExtendedModel(Alpha, Beta, Eta, Delta, C, T, S, RM, f, lb, ub, iterations, fname, varargin)

if numel(varargin) == 2
    pertubation = varargin{1};
    sampleTimes = varargin{2};
else
    pertubation = nan;
    sampleTimes = [];
end

% Draw random base omegas and costs
omega_0 = 10.^(lb + (ub-lb) * rand(size(RM)));
cost    = 1-f*rand(size(RM));

% Prepare empty lists
B = {};             % List of bacteria
age = [];           % List of ages

% Check for pertubation
if ~isnan(pertubation)
    if pertubation > 1
        P = {[], []};
    else
        P = {[], [], [], [], []};
    end
else
    P = {[]};
end

I = [];

B_samples = {};

bacteria  = nan(1, iterations);
if isnan(pertubation)
    phages = nan(1, iterations);
else
    phages = nan(numel(P), iterations);
end
diversity = nan(1, iterations);
overlap   = nan(1, iterations);
mRM       = nan(1, iterations);

% Main loop
for i = 1:iterations

    if i <= iterations

        % Draw one more species
        B = generateSpeciesExtended(B, S, RM);

        % Store the age
        age = [age; zeros(numel(B)-numel(age), 1)];

    end

    % Count the number of species
    nB = numel(B);
    nP = numel(P);

    % Set starting conditions
    if i == 1
        x0 = ones(2*nB, 1); % Everyone starts with a population of 1
    else
        b = ones(nB-sum(I), 1); % New bacteria and phages start with a population of 1
        x0 = [y(1:sum(I), end); b; abs(y(sum(I)+(1:sum(I)), end)); b];
    end

    % Check for pertubation (adding phage)
    if ~isnan(pertubation)
        if i < pertubation
            x0 = [x0; zeros(nB*(nP-1), 1)];
        elseif i == pertubation
            x0 = [x0; ones(nB*(nP-1), 1)];

            % Compute weights
            w = zeros(1, numel(RM));
            for b = 1:nB
                w(B{b}) = w(B{b}) + x0(b);
            end

            if pertubation > 1

                % Get sorting indicies
                [~, I] = sort(w, 'desc');

                % Find the three most frequent RM system
                P = {[], I(1:3)};

            else

                % Randomly give immunity to 1/5 of RM systems to each of
                % the five phages
                P = {sort(randperm(RM(end), floor(numel(RM)*0.2))), sort(randperm(RM(end), floor(numel(RM)*0.2))), sort(randperm(RM(end), floor(numel(RM)*0.2))), sort(randperm(RM(end), floor(numel(RM)*0.2))), sort(randperm(RM(end), floor(numel(RM)*0.2)))};

            end

        elseif i > pertubation
            x0 = [x0; reshape([reshape(y(2*sum(I)+1:end, end), sum(I), nP-1); repmat(b, 1, nP-1)], nB*(nP-1), 1)];
        end
    end

    % Run dynamics
    [~, ~, y] = DynamicalSystem('ExtendedModel', B, P, false, omega_0, cost, {Alpha, Beta, Eta, Delta, C, T}, x0);

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
    phages(:, i) = sum(reshape(P_end, nB, size(phages, 1)));
    D_end        = diversity(i);
    mRM(i)       = mean(cellfun(@numel, B));

    % Remove the dead species
    y(repmat(~I, nP+1, 1), :) = [];
    B(~I) = [];
    age(~I) = [];

    % Increase the age
    age = age + 1;

    overlap(i) = mean(arrayfun(@(r) sum([B{:}]==r), unique([B{:}])) - 1)/(D_end-1);

    % Sample B if required
    if any(i == sampleTimes)
        B_samples{end+1} = B;
    end

    if mod(i, 100) == 0
        save(fname, 'B', 'bacteria', 'phages', 'diversity', 'overlap', 'mRM', 'cost', 'omega_0', 'age', 'B_samples', 'C', 'Alpha', 'Beta', 'Delta', 'T', 'lb', 'ub', 'S', 'f', 'iterations', 'pertubation', 'sampleTimes')
    end

end
end