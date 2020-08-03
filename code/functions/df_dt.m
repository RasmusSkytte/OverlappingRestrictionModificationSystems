function y = df_dt(x, nB, gamma, omega, C, alpha, beta, eta, delta)

nP = (numel(x)-nB)/nB;

% Allocate y
y = zeros(size(x));

% Compute B
B = sum(x(1:nB));

% Compute eta * omega * p
eop = zeros(nB, nP);
for j = 1:nP
    eop(:, j) = eta * omega(:, :, j) * x(j*nB + (1:nB));
end

% Equations for bacteria
y(1:nB) = (gamma * (1 - B / C) - alpha - sum(eop, 2) ) .* x(1:nB);

% Phage profilation, adsorption and decay
for j = 1:nP
    y(j*nB + (1:nB)) = beta * eop(:, j) .* x(1:nB) -  x(j*nB + (1:nB)) * (eta * B + delta);
end

end