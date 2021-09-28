function [A_ij, A_RM] = computeNetwork(B)

% Compute the span of the RM systems
RM = 1:max([B{:}]);

A_ij = zeros(numel(B), numel(RM));
A_RM = zeros(numel(RM));

% Presence absence matrix
for i = 1:numel(B)
    A_ij(i, B{i}) = 1;
end

% Overlap matrix
for ri = unique([B{:}])
    for rj = unique([B{:}])
        for b = 1:numel(B)
            if ri == rj
                A_RM(ri, ri) = 1;
                continue;
            end
            if ismember(RM(ri),B{b}) && ismember(RM(rj),B{b})
                A_RM(ri, rj) = A_RM(ri, rj) + 1;
            end
        end
    end
end

if ~isempty(A_RM)
    
    % Prune RM network
    I = sum(A_RM) == 0;
    A_RM(I,:) = [];
    A_RM(:,I) = [];
    
    % Remove self links
    A_RM = A_RM - eye(size(A_RM));

    % Prune adjacency matrix network
    A_ij = A_ij(:, sum(A_ij)>0);
end
end