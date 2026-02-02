%% ALL_TOY_PROBLEMS - Main script to run and compare all toy problem models
%
% This script serves as the main entry point for running various dynamical
% systems toy problems for testing GP-SINDy and related methods.
%
% Available models:
%   1. Logistic Growth          - logistic_growth.m
%   2. Lotka-Volterra           - lotka_volterra_model.m
%   3. Reaction Kinetics        - reaction_kinetics.m
%   4. SINDy+MCMC Microglia     - sindy_mcmc_microglia.m
%   5. ESINDy+DBN Microglia     - esindy_dbn_microglia.m
%
% Author: Your Name
% Date: 2026

clear; clc; close all;

%% ========================================================================
%  SELECT WHICH MODELS TO RUN
%  ========================================================================
run_logistic = true;
run_lotka_volterra = true;
run_reaction_kinetics = true;
run_sindy_mcmc = true;
run_esindy_dbn = true;
run_comparison = true;  % Side-by-side comparison of microglia models

%% ========================================================================
%  1. LOGISTIC GROWTH
%  ========================================================================
if run_logistic
    fprintf('\n========== LOGISTIC GROWTH ==========\n');
    
    % Use default parameters
    [t_log, P_log] = logistic_growth();
    
    % Or customize:
    % params_log.r = 0.8;
    % params_log.K = 500;
    % [t_log, P_log] = logistic_growth(params_log);
end

%% ========================================================================
%  2. LOTKA-VOLTERRA PREDATOR-PREY
%  ========================================================================
if run_lotka_volterra
    fprintf('\n========== LOTKA-VOLTERRA ==========\n');
    
    % Use default parameters
    [t_lv, prey, predators] = lotka_volterra_model();
    
    % Or customize:
    % params_lv.alpha = 1.5;
    % params_lv.beta = 0.5;
    % [t_lv, prey, predators] = lotka_volterra_model(params_lv);
end

%% ========================================================================
%  3. STIFF REACTION KINETICS
%  ========================================================================
if run_reaction_kinetics
    fprintf('\n========== REACTION KINETICS ==========\n');
    
    % Use default parameters
    [t_rk, C_rk] = reaction_kinetics();
    
    % Or customize:
    % params_rk.k1 = 1e4;
    % params_rk.k2 = 1e-1;
    % [t_rk, C_rk] = reaction_kinetics(params_rk);
end

%% ========================================================================
%  4. SINDy+MCMC MICROGLIAL CELL MODEL (Amato & Arnold, 2025)
%  ========================================================================
if run_sindy_mcmc
    fprintf('\n========== SINDy+MCMC MICROGLIA ==========\n');
    
    % Use default parameters (MCMC posterior means from paper)
    [t_sindy, M1_sindy, M2_sindy, M_eq_sindy] = sindy_mcmc_microglia();
    
    % Or customize:
    % params_sindy.tspan = [0 100];  % Longer simulation
    % [t_sindy, M1_sindy, M2_sindy] = sindy_mcmc_microglia(params_sindy);
end

%% ========================================================================
%  5. ESINDy+DBN MICROGLIAL CELL MODEL (Amato & Arnold, 2025)
%  ========================================================================
if run_esindy_dbn
    fprintf('\n========== ESINDy+DBN MICROGLIA ==========\n');
    
    % Use default parameters (averaged coefficients from paper)
    [t_esindy, M1_esindy, M2_esindy] = esindy_dbn_microglia();
    
    % Or customize:
    % params_esindy.t_end = 100;  % Longer simulation
    % [t_esindy, M1_esindy, M2_esindy] = esindy_dbn_microglia(params_esindy);
end

%% ========================================================================
%  6. COMPARISON: SINDy+MCMC vs ESINDy+DBN (Figure 14 style)
%  ========================================================================
if run_comparison && run_sindy_mcmc && run_esindy_dbn
    fprintf('\n========== MODEL COMPARISON ==========\n');
    
    figure('Position', [100, 100, 1200, 500]);
    
    % SINDy+MCMC subplot
    subplot(1,2,1);
    plot(t_sindy, M1_sindy, 'k-', 'LineWidth', 2); hold on;
    plot(t_sindy, M2_sindy, 'r-', 'LineWidth', 2);
    xlabel('Time (Days)', 'FontSize', 12);
    ylabel('cells/mm^2', 'FontSize', 12);
    title('SINDy+MCMC (Continuous-Time)', 'FontSize', 14);
    legend('M1', 'M2', 'Location', 'best');
    xlim([0 50]); ylim([0 1500]);
    grid on;
    
    % ESINDy+DBN subplot
    subplot(1,2,2);
    plot(t_esindy, M1_esindy, 'k-', 'LineWidth', 2); hold on;
    plot(t_esindy, M2_esindy, 'r-', 'LineWidth', 2);
    xlabel('Time (Days)', 'FontSize', 12);
    ylabel('cells/mm^2', 'FontSize', 12);
    title('ESINDy+DBN (Discrete-Time)', 'FontSize', 14);
    legend('M1', 'M2', 'Location', 'best');
    xlim([0 50]); ylim([0 1500]);
    grid on;
    
    sgtitle('Comparison: Data-Driven Models of Microglial Cell Dynamics (Amato & Arnold 2025)', 'FontSize', 14);
    
    % Print comparison summary
    fprintf('\nModel Comparison Summary:\n');
    fprintf('  SINDy+MCMC final values:  M1 = %.2f, M2 = %.2f\n', M1_sindy(end), M2_sindy(end));
    fprintf('  ESINDy+DBN final values:  M1 = %.2f, M2 = %.2f\n', M1_esindy(end), M2_esindy(end));
    fprintf('  SINDy+MCMC equilibrium:   M1* = %.2f, M2* = %.2f\n', M_eq_sindy(1), M_eq_sindy(2));
end

%% ========================================================================
%  SUMMARY
%  ========================================================================
fprintf('\n========================================\n');
fprintf('All selected toy problems completed.\n');
fprintf('========================================\n');
