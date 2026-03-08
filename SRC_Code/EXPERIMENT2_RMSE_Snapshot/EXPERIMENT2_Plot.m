% =========================================================================
% Title: RMSE vs. Snapshots Plot for OAM-FDA MIMO Radar
% Description: 
% Generates 3 independent single-column figures.
% =========================================================================
clear; close all; clc;
load('RMSE_theta_prop.mat');
load('RMSE_phi_prop.mat'  );
load('RMSE_range_prop.mat');
load('CRLB_theta_k.mat' );
load('CRLB_phi_k.mat'   );
load('CRLB_range_k.mat' );
load('MUSIC_theta.mat'    );
load('MUSIC_phi.mat'      );
load('MUSIC_range.mat'    );



K_vec = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500];
c_prop = [0 0.4470 0.7410]; c_mus = [0.8500 0.3250 0.0980]; c_crlb = [0 0 0];
fig_size = [100, 100, 420, 360];

fig1 = figure('Name', 'RMSE_Snap_Theta', 'Position', fig_size, 'Color', 'w');
semilogy(K_vec, MUSIC_theta, '-s', 'LineWidth', 1.8, 'MarkerSize', 7, 'Color', c_mus); hold on;
semilogy(K_vec, RMSE_theta_prop, '-o', 'LineWidth', 1.8, 'MarkerSize', 7, 'Color', c_prop);
semilogy(K_vec, CRLB_theta_k, '--', 'LineWidth', 2.0, 'Color', c_crlb);
grid on; set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
xlabel('Number of Snapshots $K$', 'Interpreter', 'latex', 'FontSize', 13); ylabel('RMSE of Elevation $\theta$ (deg)', 'Interpreter', 'latex', 'FontSize', 13);
xlim([min(K_vec) max(K_vec)]); xticks(K_vec); ylim([1e-3, 1e2]);
lgd1 = legend('3-D MUSIC', 'Proposed Method', 'CRLB', 'Location', 'northeast'); set(lgd1, 'FontSize', 10, 'Interpreter', 'latex');

fig2 = figure('Name', 'RMSE_Snap_Phi', 'Position', fig_size + [450, 0, 0, 0], 'Color', 'w');
semilogy(K_vec, MUSIC_phi, '-s', 'LineWidth', 1.8, 'MarkerSize', 7, 'Color', c_mus); hold on;
semilogy(K_vec, RMSE_phi_prop, '-o', 'LineWidth', 1.8, 'MarkerSize', 7, 'Color', c_prop);
semilogy(K_vec, CRLB_phi_k, '--', 'LineWidth', 2.0, 'Color', c_crlb);
grid on; set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
xlabel('Number of Snapshots $K$', 'Interpreter', 'latex', 'FontSize', 13); ylabel('RMSE of Azimuth $\phi$ (deg)', 'Interpreter', 'latex', 'FontSize', 13);
xlim([min(K_vec) max(K_vec)]); xticks(K_vec); ylim([1e-3, 1e2]);
lgd2 = legend('3-D MUSIC', 'Proposed Method', 'CRLB', 'Location', 'northeast'); set(lgd2, 'FontSize', 10, 'Interpreter', 'latex');

fig3 = figure('Name', 'RMSE_Snap_Range', 'Position', fig_size + [900, 0, 0, 0], 'Color', 'w');
semilogy(K_vec, MUSIC_range, '-s', 'LineWidth', 1.8, 'MarkerSize', 7, 'Color', c_mus); hold on;
semilogy(K_vec, RMSE_range_prop, '-o', 'LineWidth', 1.8, 'MarkerSize', 7, 'Color', c_prop);
semilogy(K_vec, CRLB_range_k, '--', 'LineWidth', 2.0, 'Color', c_crlb);
grid on; set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
xlabel('Number of Snapshots $K$', 'Interpreter', 'latex', 'FontSize', 13); ylabel('RMSE of Range $r$ (m)', 'Interpreter', 'latex', 'FontSize', 13);
xlim([min(K_vec) max(K_vec)]); xticks(K_vec); ylim([1e-2, 1e2]);
lgd3 = legend('3-D MUSIC', 'Proposed Method', 'CRLB', 'Location', 'southwest'); set(lgd3, 'FontSize', 10, 'Interpreter', 'latex');
