function B = generateSpeciesExtended(B, S, RM)

nB = numel(B);

while numel(B) < max(S, nB+1)
    
    % If there are less than S species, add completely new species
    if numel(B) < S
        b = RM(randi(numel(RM)));
        
    else % If there are more than S species, mutate from existing
     
        % Draw from poisson distribution
        b = poissrnd(median(cellfun(@numel, B)));
        
        % Ensure we are within bounds
        b = min(numel(RM), max(1, b));
        
        % Choose random RM
        b = sort(RM(randperm(numel(RM), b)));
        
    end
    
    % Check whether it exists
    exist = false;
    for j = 1:numel(B)
        if all(ismember(b, B{j})) && all(ismember(B{j}, b))
            exist = true;
        end
    end
    if exist
        continue
    end
    
    % Store the bacterium
    B{end+1} = b;
    
end
end
