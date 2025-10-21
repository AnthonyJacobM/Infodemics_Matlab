params = struct('chi',0.048,'chi_hat',0.37,...
                'epsilon',0.33,'mu',0.1,...
                'chi_bar',0.0048,'gam',0.07,'delta',0.90,...
                'r',1.635,'s',1,'f',0);

default_dfe_ics = [ 0.175, 0.375, 0.09, 0.125, 0.4 ]; % mu = 1e-6
default_sf_ics = [ 0.125, 0.4, 0.11, 0.25, 0.1 ]; % mu = 0.25
default_lc_ics = [ 0.225, 0.475, 0.085, 0.125, 0.0 ]; % mu = 0.066

y0 = [0.9 0.05 0.02 0.01 0.02];     % initial conditions
yss = [0.85 0.06 0.03 0.02 0.04 0.01]; % example - steady state vector
t = linspace(0,10000,10001);

% Fig2F:
params.mu = 0.25;
y0 = default_sf_ics;
% Fig2C:
params.mu = 1e-6;
y0 = default_dfe_ics;
% Fig2I:
params.mu = 0.066;
y0 = default_lc_ics;
yss = [ 0.0835784742933273, 0.417846388379996, 0.13692599340057, 0.188788486675675, 0.000641150687453772 ]; % mu = 0.25
yss = [ 0.341306089623631, 0.275416966101144, 0, 0.383276944275225, 0 ]; % mu = 1e-6
yss = [ 0.241136377556403, 0.533230964751255, 0.0724586272033277, 0.00361909512384761, 0 ]; % mu = 0.066
plot_nullclines2(y0, params, t, yss, 'B');




function plot_nullclines2(y0, params, t, yss, option)
% plot_nullclines(y0, params, t, yss, option)
% -------------------------------------------------
% Plots phase portraits and nullclines for the infodemics model.
%
% INPUTS:
%   y0     : initial conditions vector
%   params : structure containing parameter values
%   t      : time vector for long-term integration
%   yss    : steady-state vector (yss = [x1_ss, x2_ss, x3_ss, x4_ss, x5_ss, ig_ss])
%   option : 'A', 'B', or 'C' to select variable pair
%
% Author: Anthony
% Date: <Date>
% -------------------------------------------------

% Extract parameter values
epsilon         = params.epsilon;
mu              = params.mu;
chi_hat         = params.chi_hat;
chi             = params.chi;
chi_bar         = params.chi_bar;
gam             = params.gam;
delta           = params.delta;
r               = params.r;
s               = params.s;
f               = params.f;

% Steady states
x1_ss = yss(1);
x2_ss = yss(2);
x3_ss = yss(3);
x4_ss = yss(4);
x5_ss = yss(5);
ig_ss = 1 - (x1_ss + x2_ss + x3_ss + x4_ss);

% Integrate system to get trajectories
[~, Y] = ode45(@(t,y) infodemics_rhs(t, y, params), t, y0);

x1_traj = Y(:,1);
x2_traj = Y(:,2);
x3_traj = Y(:,3);
x4_traj = Y(:,4);
x5_traj = Y(:,5);

n_bin = 200; % grid resolution
xlow = []; xhigh = []; ylow = []; yhigh = [];

%% CASE A: S_G vs S_B
if strcmp(option,'A') || isempty(option)
    xlab = '$S_G$';
    ylab = '$S_B$';
    xnull_lab = '$N_{S_G}$';
    ynull_lab = '$N_{S_B}$';

    x_traj = x1_traj(1:end);
    y_traj = x2_traj(1:end);

    X_SS = x1_ss; % Steady State XCoord
    Y_SS = x2_ss; % ...          YCoord

    % Axis limits
    if isempty(xlow), xlow = 0.95 * min(x_traj); end
    if isempty(xhigh), xhigh = 1.05 * max(x_traj); end
    if isempty(ylow), ylow = 0.95 * min(y_traj); end
    if isempty(yhigh), yhigh = 1.05 * max(y_traj); end

    x_array = linspace(xlow, xhigh, n_bin);
    y_array = linspace(ylow, yhigh, n_bin);

    % Nullclines
    x_null = (gam * ig_ss - x_array .* (x5_ss + (chi + mu) * x3_ss)) ./ ...
             (mu * x_array);
    y_null = (x3_ss * (chi_hat * y_array - gam)) ./ ...
             (mu * (y_array + x3_ss));

%% CASE B: S_B vs I_B
elseif strcmp(option,'B')
    xlab = '$S_B$';
    ylab = '$I_B$';
    xnull_lab = '$N_{S_B}$';
    ynull_lab = '$N_{I_B}$';

    x_traj = x2_traj(1:end);
    y_traj = x3_traj(1:end);

    X_SS = x2_ss; % Steady State XCoord
    Y_SS = x3_ss; % ...          YCoord

    if isempty(xlow), xlow = 0.95 * min(x_traj); end
    if isempty(xhigh), xhigh = 1.05 * max(x_traj); end
    if isempty(ylow), ylow = 0.95 * min(y_traj); end
    if isempty(yhigh), yhigh = 1.05 * max(y_traj); end

    x_array = linspace(xlow, xhigh, n_bin);
    y_array = linspace(ylow, yhigh, n_bin);

    % Nullclines
    x_null = (mu * x1_ss .* x_array) ./ ...
             (chi_hat * x_array - gam - mu * x1_ss);
    y_null = (gam + epsilon * (x1_ss + ig_ss) - x5_ss * x1_ss ./ y_array) / chi_hat;
    %y_null = (gam * y_array + epsilon * y_array * (x1_ss + ig_ss) - chi_bar * y_array * x4_ss) ./ (y_array * (1 - delta) * chi);
    
%% CASE C: I_B vs V
elseif strcmp(option,'C')
    xlab = '$I_B$';
    ylab = '$V$';
    xnull_lab = '$N_{I_B}$';
    ynull_lab = '$N_{V}$';

    x_traj = x3_traj(2000:end);
    y_traj = x4_traj(2000:end);

    X_SS = x3_ss; % Steady State XCoord 
    Y_SS = x4_ss; % ...          YCoord

    if isempty(xlow), xlow = 0.95 * min(x_traj); end
    if isempty(xhigh), xhigh = 1.05 * max(x_traj); end
    if isempty(ylow), ylow = 0.95 * min(y_traj); end
    if isempty(yhigh), yhigh = 1.05 * max(y_traj); end

    x_array = linspace(xlow, xhigh, n_bin);
    y_array = linspace(ylow, yhigh, n_bin);

    % Nullclines
    x_null = ones(1, n_bin) * ...
        (gam + epsilon * (x1_ss + ig_ss) - chi_hat * x2_ss) / ...
        ((1 - delta) * chi); % V = rhs
    y_null = (x5_ss * x1_ss) ./ ((1 - delta) * chi * x_array); % V = rhs

end

%% --- Plot results ---
figure; hold on; box on;
plot(x_traj, y_traj, 'b', 'LineWidth', 1.2);
plot(x_array, x_null, 'r--', 'LineWidth', 1.5, 'DisplayName', xnull_lab);
plot(y_null, y_array, 'k-.', 'LineWidth', 1.5, 'DisplayName', ynull_lab);
%scatter(X_SS, Y_SS, 70, 'ko', 'filled', 'DisplayName', 'Steady State');
xlabel(xlab, 'Interpreter', 'latex');
ylabel(ylab, 'Interpreter', 'latex');
legend('Interpreter','latex','Location','best');
title([strcat('Phase Portrait with Nullclines (', xlab, ' vs ', ylab, ')')]);
grid on;

end

function dydt = infodemics_rhs(~,y,params)
    sg  = y(1);
    sb  = y(2);
    ib  = y(3);
    v   = y(4);
    phi = y(5);

    gam = params.gam; chi = params.chi; chi_bar = params.chi_bar; chi_hat = params.chi_hat;
    epsilon = params.epsilon; mu = params.mu; r = params.r; s = params.s;

    ig = 1 - (sg + sb + ib + v);

    sgp  = gam*ig - phi*sg - chi*ib*sg - mu*sg*(sb + ib);
    sbp  = gam*ib + mu*sg*(sb + ib) - chi_hat*ib*sb;
    ibp  = ib*(-gam + chi_hat*sb + chi_bar*v - epsilon*(sg + ig));
    vp   = phi*sg - chi_bar*ib*v;
    phip = s*phi*(1 - phi)*(ig + ib - r*v);

    dydt = [sgp; sbp; ibp; vp; phip];

end
