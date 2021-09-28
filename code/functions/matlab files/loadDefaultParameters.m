% Define default parameters
C       = 1e8;  % Carrying capacity
Eta     = 1e-8; % Adsorption rate
Alpha   = 0.2;  % Dilution rate of bacteria
Beta    = 100;  % Burst size of phages
Delta   = 0.2;  % Decay rate of phages
T       = 1e3;  % Time between addition of species
lb      = -4;   % Lower bound for omega
ub      = 0;    % Upper bound for omega
N       = 1e4;  % Number of points in distribution clouds
S       = 5;    % Min number of species

% Compute cost factor
f = 0.1;

iterations = 1e6; 	% Number of iterations