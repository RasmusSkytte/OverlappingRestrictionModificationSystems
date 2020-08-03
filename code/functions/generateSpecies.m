function [gamma, omega, n] = generateSpecies(nRM, lb, ub, f, S, N)

% Allocate
gamma = nan(N, 1);
omega = nan(N, 1);
n     = nan(N, 1);

nB = numel(nRM);

k = 1;
while numel(nRM) < max(S, nB+1)
    
    % If there are less than S species, add completely new species
    if numel(nRM) < S
        m = 1;
        
    else % If there are more than S species, mutate from existing
        
        % Draw from poisson distribution
        m = poissrnd(median(nRM));
       
        % Ensure we are within bounds
        m = max(1, m);
        
    end
    
    % Sample the distribution
    gamma(k) = prod(1-f*rand(m, 1));
    omega(k) = prod(10.^(lb + (ub-lb) * rand(m, 1)));
    n(k)     = m;
    nRM = [nRM; m];
    
    k = k + 1;
    
end

end
