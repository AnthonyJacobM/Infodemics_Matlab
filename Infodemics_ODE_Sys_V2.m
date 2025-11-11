%% =================== MAIN SCRIPT: multiCaseDynamics.m ===================
clc; clear; close all;

% === Load Default Parameters ===
params_default = struct('gam',0.07,'chi',0.048,'chi_bar',0.0048,'chi_hat',0.37,...
                        'epsilon',0.33,'mu',0.066,'r',1.635,'s',1);

% === Default Steady States and ICs ===
default_dfe_ics = [0.175, 0.375, 0.09, 0.125, 0.4];
default_sf_ics  = [0.125, 0.4, 0.11, 0.25, 0.1];
default_lc_ics  = [0.225, 0.475, 0.085, 0.125, 0.0];
default_lc_ss   = [0.213899947767176, 0.486438741030547, 0.0716742222372843, 0.103285708601247, 2.1966e-17];
default_r_bp_ss = [0.000571591255633368, 0.189492234608488, 0.197046892184315, 0.604330828414679, 1];

% === Select Case (A–H) ===
caseChoice = 'H';  % Change this to run a different case ('A'...'H')

% === Get Parameters and Initial Conditions ===
[params, y0, figNames] = getCaseParams(caseChoice, params_default, ...
    default_dfe_ics, default_sf_ics, default_lc_ics, default_lc_ss, default_r_bp_ss);

% === Solve ODE System ===
tspan = [0 15000];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t, y] = ode45(@(t,y) systemODE(t,y,params), tspan, y0, opts);

% === Compute Derived Quantities ===
[ig, I, G, B, Re] = derivedQuantities(y, params);

% === Plotting (Unified Style) ===
plotStateVariables(t, y, figNames{1});
plotRe(t, Re, figNames{2});
plotPhasePlane(y, figNames{3}, caseChoice);

disp(['✅ Case ' caseChoice ' completed successfully.']);

%% ================= Supporting Functions =================

function [params, y0, figNames] = getCaseParams(caseChoice, params_default, dfe, sf, lc, lc_ss, r_bp_ss)
    switch upper(caseChoice)
        case 'A'
            params = params_default; params.mu = 1e-6;
            y0 = [0.15, 0.375, 0.08, 0.12, 0.4];  % y0_2A
            figNames = {'Fig2A','Fig2B','Fig2C'};
        case 'B'
            params = params_default; params.mu = 0.066;
            y0 = lc;
            figNames = {'Fig2D','Fig2E','Fig2F'};
        case 'C'
            params = params_default; params.mu = 0.25;
            y0 = sf;
            figNames = {'Fig2G','Fig2H','Fig2I'};
        case 'D'
            params = params_default; params.mu = 0.10; params.r = 0.34022;
            y0 = r_bp_ss - rand(1,5)*1e-2;
            figNames = {'Fig3A','Fig3B','Fig3C'};
        case 'E'
            params = params_default; params.mu = 0.066; params.epsilon = 0.4438;
            y0 = r_bp_ss - rand(1,5)*1e-4;
            figNames = {'Fig3D','Fig3E','Fig3F'};
        case 'F'
            params = params_default; params.r = 0.0908; params.mu = 0.066; params.epsilon = 0.9033;
            y0 = r_bp_ss - rand(1,5)*1e-4;
            figNames = {'Fig3G','Fig3H','Fig3I'};
        case 'G'
            params = params_default; params.r = 1.635; params.mu = 0.0660611979693606; params.epsilon = 0.33;
            yss = [0.209937863406593, 0.482210866198904, 0.0707320122282017, 0.116818487531245, 0.000188920052159354];
            y0 = yss + rand(1,5)*0.01;
            figNames = {'Fig4C','Fig4D','Fig4E'}; 
        case 'H'
            params = params_default; params.r = 0.0908; params.mu = 0.066; params.epsilon = 0.9033;
            y0 = [0.000571591255633368, 0.189492234608488, 0.197046892184315, 0.604330828414679, 0.1*rand(1)];
            figNames = {'Fig6C','Fig6D','Fig6E'};
        otherwise
            error('Invalid case selection. Choose A–H.');
    end
end

function dydt = systemODE(~,y,p)
    sg = y(1); sb = y(2); ib = y(3); v = y(4); phi = y(5);
    ig = 1 - (sg + sb + ib + v);
    dydt = [
        p.gam*ig - phi*sg - p.chi*ib*sg - p.mu*sg*(sb + ib);
        p.gam*ib + p.mu*sg*(sb + ib) - p.chi_hat*ib*sb;
        ib*(-p.gam + p.chi_hat*sb + p.chi_bar*v - p.epsilon*(sg + ig));
        phi*sg - p.chi_bar*ib*v;
        p.s*phi*(1 - phi)*(ig + ib - p.r*v)
    ];
end

function [ig, I, G, B, Re] = derivedQuantities(y,p)
    sg = y(:,1); sb = y(:,2); ib = y(:,3); v = y(:,4);
    ig = 1 - (sg + sb + ib + v);
    I = ig + ib; G = sg + ig; B = sb + ib;
    Re = (p.chi_hat.*sb + p.chi_bar.*v) ./ (p.gam + p.epsilon.*(sg + ig));
end

function plotStateVariables(t, y, figName)
    figure('Name',figName,'Position',[100 100 900 600]);
    plot(t, y(:,1),'b','LineWidth',1.5); hold on;
    plot(t, y(:,2),'r','LineWidth',1.5);
    plot(t, y(:,3),'r--','LineWidth',1.5);
    plot(t, y(:,4),'Color',[0 0.5 0],'LineWidth',1.5);
    plot(t, y(:,5),'k','LineWidth',1.5);
    xlabel('Time'); ylabel('State Variables'); ylim([0 1]); grid on;
    legend('S_G','S_B','I_B','V','\phi','Location','best');
end

function plotRe(t, Re, figName)
    figure('Name',figName,'Position',[100 100 900 600]);
    plot(t, Re, 'b','LineWidth',1.5); xlabel('Time'); ylabel('R_e(t)');
    ylim([0.75,1.25]); grid on;
end

function plotPhasePlane(y, figName, caseChoice)
    figure('Name',figName,'Position',[100 100 900 600]);
    switch upper(caseChoice)
        case {'A','B','E','H'}
            x = y(:,2); yvar = y(:,3); xlabel('S_B'); ylabel('I_B');
        case 'C'
            x = y(:,1); yvar = y(:,2); xlabel('S_G'); ylabel('S_B');
        case 'D'
            x = y(:,3); yvar = y(:,4); xlabel('I_B'); ylabel('V');
        case 'F'
            x = y(:,1); yvar = y(:,2); xlabel('S_G'); ylabel('S_B');
        case 'G'
            x = y(:,2); yvar = y(:,3); xlabel('S_B'); ylabel('I_B');
        otherwise
            x = y(:,3); yvar = y(:,4); xlabel('I_B'); ylabel('V');
    end
    plot(x, yvar, 'LineWidth', 1.5, 'Color', 'b');
    grid on; title(['Phase Trajectory: ' figName]);
end
