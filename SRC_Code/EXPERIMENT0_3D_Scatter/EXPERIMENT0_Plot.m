% =========================================================================
% Title: 3-D Scatter Plot for OAM-FDA MIMO Radar
% Description: 
% 1. Visualizes 3-D localization and automatic pairing.
% 2. Runs the actual CP-ALS algorithm for authentic visualization.
% =========================================================================

clear; close all; clc;
rng(2026); % Set random seed for reproducibility

%% 1. Radar & Target Parameters
c = 3e8; fc = 10e9; lambda = c / fc;
N_rx = 12; d_rx = lambda / 2;
l_modes = [-3, -2, -1, 0, 1, 2, 3]'; P_oam = length(l_modes); R_uca = 2 * lambda;
M_fda = 10; delta_f = 20e3;
K_snap = 500; SNR_dB = 15;

theta_true = [9.5, 21.3, 32.7, 46.1, 63.3, 74.8] * pi / 180;
phi_true   = [-21.2,-33.6,-42.7, 18.8, -21.3, 51.4] * pi / 180;
range_true = [201, 1903, 822, 2251, 4869, 6507];
L = length(range_true);

%% 2. Generate True Manifolds 
k_wave = 2 * pi / lambda;
A_true = zeros(N_rx, L); B_true = zeros(P_oam, L); C_true = zeros(M_fda, L);
for i = 1:L
    A_true(:,i) = exp(1j * 2 * pi * d_rx * (0:N_rx-1)' * cos(theta_true(i)) / lambda);
    B_true(:,i) = besselj(l_modes, k_wave * R_uca * sin(theta_true(i))) .* exp(1j * l_modes * phi_true(i));
    C_true(:,i) = exp(-1j * 4 * pi * (0:M_fda-1)' * delta_f * range_true(i) / c);
end
U_true = khatri_rao(C_true, khatri_rao(B_true, A_true));
signal_pwr_base = sum(sum(abs(U_true).^2)) / size(U_true,1); % Exact power alignment
noise_var_actual = signal_pwr_base / (10^(SNR_dB / 10));

%% 3. Generate Tensor & Run CP-ALS
S_true = (randn(K_snap, L) + 1j * randn(K_snap, L)) / sqrt(2);
Y_clean = zeros(N_rx, P_oam, M_fda, K_snap);
for i = 1:L
    t4d = reshape(A_true(:, i) * B_true(:, i).', [N_rx, P_oam, 1, 1]);
    t4d = reshape(t4d(:) * C_true(:, i).', [N_rx, P_oam, M_fda, 1]);
    t4d = reshape(t4d(:) * S_true(:, i).', [N_rx, P_oam, M_fda, K_snap]);
    Y_clean = Y_clean + t4d;
end

noise = (randn(size(Y_clean)) + 1j * randn(size(Y_clean))) / sqrt(2);
Y_noisy = Y_clean + sqrt(noise_var_actual) * noise;

disp('Running authentic CP-ALS for 3D Scatter Plot...');
[A_est, B_est, C_est, ~] = cp_als_4d_custom(Y_noisy, L, 100, 1e-5, A_true, B_true, C_true, S_true);

theta_est = zeros(1, L); phi_est = zeros(1, L); range_est = zeros(1, L);
for i = 1:L
    a_v = A_est(:, i);
    ph_z = angle((a_v(1:end-1)' * a_v(1:end-1)) \ (a_v(1:end-1)' * a_v(2:end)));
    theta_est(i) = acos(max(min(ph_z * lambda / (2 * pi * d_rx), 1), -1));
    
    b_v = B_est(:, i);
    W_l = besselj(l_modes, k_wave * R_uca * sin(theta_est(i)));
    b_comp = b_v .* sign(W_l); 
    phi_est(i) = angle((b_comp(1:end-1)' * b_comp(1:end-1)) \ (b_comp(1:end-1)' * b_comp(2:end)));
    
    c_v = C_est(:, i);
    meas_phase = angle((c_v(1:end-1)' * c_v(1:end-1)) \ (c_v(1:end-1)' * c_v(2:end)));
    if meas_phase > 0, true_phase = meas_phase - 2*pi; else, true_phase = meas_phase; end
    range_est(i) = true_phase * c / (-4 * pi * delta_f);
end
[tE, pE, rE] = alignParameters(theta_est*180/pi, phi_est*180/pi, range_est, theta_true*180/pi, phi_true*180/pi, range_true);

%% 4. Plotting with 3-Plane Drop Lines
fig = figure('Name', '3D_Scatter_MultiPlane', 'Position', [100, 100, 700, 550], 'Color', 'w');
hold on; grid on; box on;

theta_t_deg = theta_true * 180/pi; 
phi_t_deg = phi_true * 180/pi;

% Define the physical boundaries (the "walls") for projection
x_min = 0;      % Range wall
y_min = -60;    % Azimuth wall (slightly wider than -40)
z_min = 0;      % Elevation floor

for i = 1:L
    % 1. Drop-line to XY plane (Bottom Floor: Z = z_min)
    plot3([range_true(i), range_true(i)], [phi_t_deg(i), phi_t_deg(i)], [z_min, theta_t_deg(i)], ...
        'Color', [0.8 0.8 0.8], 'LineStyle', '--', 'LineWidth', 1.0, 'HandleVisibility','off');
    
    % 2. Drop-line to XZ plane (Back Wall: Y = y_min)
    plot3([range_true(i), range_true(i)], [y_min, phi_t_deg(i)], [theta_t_deg(i), theta_t_deg(i)], ...
        'Color', [0.8 0.8 0.8], 'LineStyle', '--', 'LineWidth', 1.0, 'HandleVisibility','off');
        
    % 3. Drop-line to YZ plane (Side Wall: X = x_min)
    plot3([x_min, range_true(i)], [phi_t_deg(i), phi_t_deg(i)], [theta_t_deg(i), theta_t_deg(i)], ...
        'Color', [0.8 0.8 0.8], 'LineStyle', '--', 'LineWidth', 1.0, 'HandleVisibility','off');
end

h1 = scatter3(range_true, phi_t_deg, theta_t_deg, 160, 'MarkerEdgeColor', [0.8500 0.3250 0.0980], 'MarkerFaceColor', 'none', 'LineWidth', 2.0);
h2 = scatter3(rE, pE, tE, 120, 'MarkerEdgeColor', [0 0.4470 0.7410], 'Marker', 'hexagram', 'LineWidth', 2.0);


set(gca, 'XLim', [x_min, 7500], 'YLim', [y_min, 60], 'ZLim', [z_min, 90]);
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman', 'LineWidth', 1.2);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);

xlabel('Target Range $r$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Azimuth $\phi$ (deg)', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('Elevation $\theta$ (deg)', 'Interpreter', 'latex', 'FontSize', 14);

lgd = legend([h1, h2], 'True Targets', 'Estimated Targets', 'Location', 'northeast');
set(lgd, 'FontSize', 12, 'Interpreter', 'latex', 'EdgeColor', [0.5 0.5 0.5]);
view(35, 25);
