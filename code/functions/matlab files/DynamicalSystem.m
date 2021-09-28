function [B_end, P_end, y, t] = DynamicalSystem(model, B, P, varargin)

originalModel = false;
extendedModel = false;
if strcmpi(model, 'OriginalModel')
    originalModel = true;
elseif strcmpi(model, 'ExtendedModel')
    extendedModel = true;
end

userDefinedOmegas = false;
if numel(varargin) >= 2
    if isvector(varargin{2})
        omega_0 = varargin{2};
        userDefinedOmegas = true;
    end
end

userDefinedCost = false;
if numel(varargin) >= 3
    if isvector(varargin{3})
        cost = varargin{3};
        userDefinedCost = true;
    end
end

userDefinedParams = false;
if numel(varargin) >= 4
    if iscell(varargin{4})
        alpha   = varargin{4}{1};
        beta    = varargin{4}{2};
        eta     = varargin{4}{3};
        delta   = varargin{4}{4};
        C       = varargin{4}{5};
        T_end   = varargin{4}{6};
        userDefinedParams = true;
    end
end

userDefinedStartingPopulations = false;
if numel(varargin) >= 5
    if isvector(varargin{5})
        x0 = varargin{5};
        userDefinedStartingPopulations = true;
    end
end

disableNormControl = false;
if numel(varargin) >= 6
    if varargin{6}
       disableNormControl = true;
    end
end

% Define variables
if ~userDefinedParams
    C      = 1e8;
    alpha  = 0.2;
    beta   = 25;
    eta    = 1/C;
    delta  = 0.2;
    T_end  = 1e4;
    RM     = unique([B{:}]);
    RM     = 1:RM(end);
end

% Draw random omegas for the RMs
if ~userDefinedOmegas
    lb = -4;
    ub = 0;
    omega_0 = 10.^(lb + (ub-lb) * rand(size(RM)));
end

% Set the cost of the RMs
if ~userDefinedCost
    cost = 1 - rand(size(RM));
end

% Get the omegas and gammas to use in the
if originalModel
    [gamma, omega] = computeGammaOmega(B, P, cost, omega_0, 'Original');
    nB = numel(cost);
    nP = max(1, numel(P));
elseif extendedModel
    [gamma, omega] = computeGammaOmega(B, P, cost, omega_0, 'Extended');
    nB = numel(B);
    nP = numel(P);
end


% Define the dynamical populations
if ~userDefinedStartingPopulations
    b = ones(nB, 1)/nB * delta/(eta*beta); % Each bacterium has equal amounts of biomass initially
    p = ones(nB*nP, 1)/(nB*nP) * mean(gamma)/eta; % Each phage has equal amounts of biomass initially
    x0 = [b; p];
end

fun = @(t, x)df_dt(x, nB, gamma, omega, C, alpha, beta, eta, delta);

options = odeset('NonNegative', find(x0 > 0), 'Normcontrol', 'on');
if disableNormControl
    options.NormControl = 'off';
end

% Run the ode solver
sol = ode45(fun, [0 T_end], x0, options);
t = sol.x;
y = sol.y;

B_end = y(1:nB, end);
P_end = sum(reshape(y((nB+1):end, end), nB, nP))';

end