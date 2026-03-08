% =========================================================================
% Title: RMSE vs. SNR Plot for OAM-FDA MIMO Radar
% Description: 
% Generates 3 independent single-column figures.
% =========================================================================
clear; close all; clc;
load('RMSE_theta_prop.mat');
load('RMSE_phi_prop.mat'  );
load('RMSE_range_prop.mat');
load('CRLB_theta_snr.mat' );
load('CRLB_phi_snr.mat'   );
load('CRLB_range_snr.mat' );
load('MUSIC_theta.mat'    );
load('MUSIC_phi.mat'      );
load('MUSIC_range.mat'    );



SNR_dB = -10:5:30; 
c_prop = [0 0.4470 0.7410]; c_mus = [0.8500 0.3250 0.0980]; c_crlb = [0 0 0];
fig_size = [100, 100, 420, 360];

fig1 = figure('Name', 'RMSE_SNR_Theta', 'Position', fig_size, 'Color', 'w');
semilogy(SNR_dB, MUSIC_theta, '-s', 'LineWidth', 1.8, 'MarkerSize', 7, 'Color', c_mus); hold on;
semilogy(SNR_dB, RMSE_theta_prop, '-o', 'LineWidth', 1.8, 'MarkerSize', 7, 'Color', c_prop);
semilogy(SNR_dB, CRLB_theta_snr, '--', 'LineWidth', 2.0, 'Color', c_crlb);
grid on; set(gca, 'YScale', 'log', 'FontSize', 11, 'FontName', 'Times New Roman');
xlabel('SNR (dB)', 'Interpreter', 'latex', 'FontSize', 13); ylabel('RMSE of Elevation $\theta$ (deg)', 'Interpreter', 'latex', 'FontSize', 13);
xlim([min(SNR_dB) max(SNR_dB)]); ylim([1e-3, 1e2]);
lgd1 = legend('3-D MUSIC', 'Proposed Method', 'CRLB', 'Location', 'northeast'); set(lgd1, 'FontSize', 10, 'Interpreter', 'latex');

fig2 = figure('Name', 'RMSE_SNR_Phi', 'Position', fig_size + [450, 0, 0, 0], 'Color', 'w');
semilogy(SNR_dB, MUSIC_phi, '-s', 'LineWidth', 1.8, 'MarkerSize', 7, 'Color', c_mus); hold on;
semilogy(SNR_dB, RMSE_phi_prop, '-o', 'LineWidth', 1.8, 'MarkerSize', 7, 'Color', c_prop);
semilogy(SNR_dB, CRLB_phi_snr, '--', 'LineWidth', 2.0, 'Color', c_crlb);
grid on; set(gca, 'YScale', 'log', 'FontSize', 11, 'FontName', 'Times New Roman');
xlabel('SNR (dB)', 'Interpreter', 'latex', 'FontSize', 13); ylabel('RMSE of Azimuth $\phi$ (deg)', 'Interpreter', 'latex', 'FontSize', 13);
xlim([min(SNR_dB) max(SNR_dB)]); ylim([1e-3, 1e2]);
lgd2 = legend('3-D MUSIC', 'Proposed Method', 'CRLB', 'Location', 'northeast'); set(lgd2, 'FontSize', 10, 'Interpreter', 'latex');

fig3 = figure('Name', 'RMSE_SNR_Range', 'Position', fig_size + [900, 0, 0, 0], 'Color', 'w');
semilogy(SNR_dB, MUSIC_range, '-s', 'LineWidth', 1.8, 'MarkerSize', 7, 'Color', c_mus); hold on;
semilogy(SNR_dB, RMSE_range_prop, '-o', 'LineWidth', 1.8, 'MarkerSize', 7, 'Color', c_prop);
semilogy(SNR_dB, CRLB_range_snr, '--', 'LineWidth', 2.0, 'Color', c_crlb);
grid on; set(gca, 'YScale', 'log', 'FontSize', 11, 'FontName', 'Times New Roman');
xlabel('SNR (dB)', 'Interpreter', 'latex', 'FontSize', 13); ylabel('RMSE of Range $r$ (m)', 'Interpreter', 'latex', 'FontSize', 13);
xlim([min(SNR_dB) max(SNR_dB)]); ylim([1e-2, 1e2]);
lgd3 = legend('3-D MUSIC', 'Proposed Method', 'CRLB', 'Location', 'northeast'); set(lgd3, 'FontSize', 10, 'Interpreter', 'latex');

