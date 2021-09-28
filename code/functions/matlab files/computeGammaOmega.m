function [gamma, omega] = computeGammaOmega(B, P, cost, omega_0, model)

if strcmpi(model,'Original')
    nB = numel(cost);
    
    gamma = cost;
    
    omega = repmat(omega_0,1,nB);
    omega = omega - diag(diag(omega)) + eye(nB);
    
elseif strcmpi(model,'Extended')
    nB = numel(B);
    nP = max(1,numel(P));
    
    % Get a list of the current RM systems
    RM = unique([B{:}]);
    if ~isempty(RM)
        RM = 1:RM(end);
    end
    
    % Compute the growth rates
    gamma = nan(nB,1);
    for i = 1:nB
        gamma(i) = prod(cost(B{i}));
    end
    
    % Generate a matrix of how the phages are methylated
    meth = cell(nB,nP);
    %           ^ methylated by bacterium
    %               ^ which phage is it?
    for i = 1:nB
        for j = 1:nP
            meth{i,j} = B{i}(~ismember(B{i},P{j}));
        end
    end
    
    % Generate a matrix of functional RM systems based on which bacteria
    % the phage has been methylated by.
    RM_eff = cell(nB, nB, nP);
    %             ^ this bacteriums percieved defence
    %                 ^ phage has been methylated by
    %                      ^ which phage is it?
    for i = 1:nB
        % Methylated phages
        for j = 1:nB
            for k = 1:nP
                RM_eff{i,j,k} = B{i}(~ismember(B{i},union(B{j},P{k})));
            end
        end
    end
    
    % Compute the RM strength against the phage phenotypes
    omega = nan(nB, nB, nP);
    %           ^ this bacteriums percieved defence
    %               ^ phage has been methylated by
    %                   ^ which phage is it?
    for i = 1:nB
        for j = 1:nB
            for k = 1:nP
                omega(i,j,k) = prod(omega_0(RM_eff{i,j,k}));
            end
        end
    end
    
else
    error('Model not known')
end
end