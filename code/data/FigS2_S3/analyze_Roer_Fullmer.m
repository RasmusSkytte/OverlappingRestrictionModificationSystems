% Analyze data

% The Patchy Distribution of Restriction�Modification System Genes and the Conservation of Orphan Methyltransferases in Halobacteria
im = importdata('Roer_Fullmer/Panel1.png');

% Define grid coordinates
x0 = 4;
xN = 333;
x  = round(linspace(x0, xN, 49));

y0 = 5;
yN = 785;
y  = round(linspace(y0, yN, 120));

% Define grid
[X, Y] = meshgrid(x, y);

% Sample the figure
A_ij = arrayfun(@(x,y)mean(im.cdata(y, x, :), 3) < 100, X, Y);

% Load panel 2
im = importdata('Roer_Fullmer/Panel2.png');

% Define grid coordinates
x0 = 4;
xN = 388;
x  = round(linspace(x0, xN, 49));

y0 = 5;
yN = 726;
y  = round(linspace(y0, yN, 97));

% Define grid
[X, Y] = meshgrid(x, y);

% Sample the figure
A_ij = [A_ij; arrayfun(@(x,y)mean(im.cdata(y, x, :),3) < 100, X, Y)];

m   = 48; % Identified as RM systems (26 for conservative estimate, 48 for full data set)
A_ij_1 = A_ij(:, 1:m);

% Is the Evolution of Salmonella entericasubsp.entericaLinked to Restriction-Modification Systems?
im = importdata('Roer_Fullmer/Panel3.tif');

% Define grid coordinates
x0 = 1138;
xN = 3766;
x  = round(linspace(x0, xN, 113));

y0 = 27;
yN = 3352;
y  = round(linspace(y0, yN, 221));

% Define grid
[X, Y] = meshgrid(x, y);

% Sample the figure
A_ij_2 = arrayfun(@(x,y)im(y, x, 1) < 100, X, Y);

% Combine into a single matrix
A_ij = [A_ij_1 zeros(size(A_ij_1, 1), size(A_ij_2, 2)); zeros(size(A_ij_2, 1), size(A_ij_1, 2)) A_ij_2];

% Label each subset
s = [ones(size(A_ij_1, 1), 1) * 1;
     ones(size(A_ij_2, 1), 1) * 2];

% Save data
save('data_Roer_Fullmer.mat', 'A_ij', 's');
