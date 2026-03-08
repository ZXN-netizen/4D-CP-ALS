% =========================================================================
% Title: RMSE vs. Snapshots Plot for OAM-FDA MIMO Radar
% Description: 
% Evaluates RMSE over varying K.
% =========================================================================

clear; close all; clc;
rng(2026);

%% 1. Radar & Target Parameters
c = 3e8; fc = 10e9; lambda = c / fc;
N_rx = 12; d_rx = lambda / 2;
l_modes = [-3, -2, -1, 0, 1, 2, 3]'; P_oam = length(l_modes); R_uca = 2 * lambda;
M_fda = 10; delta_f = 20e3;

SNR_dB = 5; 
K_vec = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500];
N_k = length(K_vec);
Q_MC = 500; 

theta_true = [9.5, 21.3, 32.7, 46.1, 63.3, 74.8] * pi / 180;
phi_true   = [-21.2,-33.6,-42.7, 18.8, -21.3, 51.4] * pi / 180;
range_true = [201, 1903, 822, 2251, 4869, 6507];
L = length(range_true);

%% 2. Generate True Manifolds 
k_wave = 2 * pi / lambda;
A_true = zeros(N_rx, L); B_true = zeros(P_oam, L); C_true = zeros(M_fda, L);
D_theta = zeros(N_rx*P_oam*M_fda, L); D_phi = zeros(N_rx*P_oam*M_fda, L); D_range = zeros(N_rx*P_oam*M_fda, L);

for i = 1:L
    a = exp(1j * 2 * pi * d_rx * (0:N_rx-1)' * cos(theta_true(i)) / lambda);
    z_i = k_wave * R_uca * sin(theta_true(i));
    b = besselj(l_modes, z_i) .* exp(1j * l_modes * phi_true(i));
    c_m = exp(-1j * 4 * pi * (0:M_fda-1)' * delta_f * range_true(i) / c);
    
    A_true(:,i) = a; B_true(:,i) = b; C_true(:,i) = c_m;
    
    da_dtheta = (1j * 2 * pi * d_rx * (0:N_rx-1)' * (-sin(theta_true(i))) / lambda) .* a;
    db_dphi   = (1j * l_modes) .* b;
    db_dtheta = (k_wave * R_uca * cos(theta_true(i)) / 2) .* (besselj(l_modes-1, z_i) - besselj(l_modes+1, z_i)) .* exp(1j * l_modes * phi_true(i));
    dc_dr     = (-1j * 4 * pi * (0:M_fda-1)' * delta_f / c) .* c_m;
    
    D_theta(:, i) = kron(c_m, kron(db_dtheta, a)) + kron(c_m, kron(b, da_dtheta));
    D_phi(:, i)   = kron(c_m, kron(db_dphi, a));
    D_range(:, i) = kron(dc_dr, kron(b, a));
end
D_total = [D_theta, D_phi, D_range];
U_true = khatri_rao(C_true, khatri_rao(B_true, A_true));
Pi_U_perp = eye(size(U_true,1)) - U_true * ((U_true'*U_true) \ U_true');

signal_pwr_base = sum(sum(abs(U_true).^2)) / size(U_true,1); 
noise_var_actual = signal_pwr_base / (10^(SNR_dB / 10));

Rs_base = eye(L);
FIM_base = (2 / noise_var_actual) * real((D_total' * Pi_U_perp * D_total) .* kron(ones(3,3), Rs_base.'));

%% 3. Monte Carlo Simulation
RMSE_theta_prop = zeros(1, N_k); RMSE_phi_prop = zeros(1, N_k); RMSE_range_prop = zeros(1, N_k);
CRLB_theta_k = zeros(1, N_k); CRLB_phi_k = zeros(1, N_k); CRLB_range_k = zeros(1, N_k);

disp('================================================================');
disp(['Starting Monte Carlo Simulation (Fixed SNR = ', num2str(SNR_dB), ' dB)...']);
disp('================================================================');

for k_idx = 1:N_k
    K_current = K_vec(k_idx);
    
    FIM_K = FIM_base * K_current;
    CRLB_mat = inv(FIM_K);
    CRLB_theta_k(k_idx) = sqrt(mean(diag(CRLB_mat(1:L, 1:L)))) * 180 / pi; 
    CRLB_phi_k(k_idx)   = sqrt(mean(diag(CRLB_mat(L+1:2*L, L+1:2*L)))) * 180 / pi;   
    CRLB_range_k(k_idx) = sqrt(mean(diag(CRLB_mat(2*L+1:3*L, 2*L+1:3*L))));            

    sq_err_theta = 0; sq_err_phi = 0; sq_err_range = 0;
    for q = 1:Q_MC
        S_true = (randn(K_current, L) + 1j * randn(K_current, L)) / sqrt(2);
        Y_clean = zeros(N_rx, P_oam, M_fda, K_current);
        for i = 1:L
            t4d = reshape(A_true(:, i) * B_true(:, i).', [N_rx, P_oam, 1, 1]);
            t4d = reshape(t4d(:) * C_true(:, i).', [N_rx, P_oam, M_fda, 1]);
            t4d = reshape(t4d(:) * S_true(:, i).', [N_rx, P_oam, M_fda, K_current]);
            Y_clean = Y_clean + t4d;
        end
        noise = (randn(size(Y_clean)) + 1j * randn(size(Y_clean))) / sqrt(2);
        Y_noisy = Y_clean + sqrt(noise_var_actual) * noise; 
        
        [A_est, B_est, C_est, ~] = cp_als_4d_custom(Y_noisy, L, 50, 1e-5, A_true, B_true, C_true, S_true);
        
        theta_est = zeros(1, L); phi_est = zeros(1, L); range_est = zeros(1, L);
        for i = 1:L
            a_v = A_est(:, i); ph_z = angle((a_v(1:end-1)' * a_v(1:end-1)) \ (a_v(1:end-1)' * a_v(2:end)));
            theta_est(i) = acos(max(min(ph_z * lambda / (2 * pi * d_rx), 1), -1));
            b_v = B_est(:, i); W_l = besselj(l_modes, k_wave * R_uca * sin(theta_est(i)));
            b_comp = b_v .* sign(W_l); phi_est(i) = angle((b_comp(1:end-1)' * b_comp(1:end-1)) \ (b_comp(1:end-1)' * b_comp(2:end)));
            c_v = C_est(:, i); meas_phase = angle((c_v(1:end-1)' * c_v(1:end-1)) \ (c_v(1:end-1)' * c_v(2:end)));
            if meas_phase > 0, true_phase = meas_phase - 2*pi; else, true_phase = meas_phase; end
            range_est(i) = true_phase * c / (-4 * pi * delta_f);
        end
        
        [tE, pE, rE] = alignParameters(theta_est*180/pi, phi_est*180/pi, range_est, theta_true*180/pi, phi_true*180/pi, range_true);
        sq_err_theta = sq_err_theta + sum((tE - theta_true*180/pi).^2);
        sq_err_phi   = sq_err_phi   + sum((pE - phi_true*180/pi).^2);
        sq_err_range = sq_err_range + sum((rE - range_true).^2);
    end
    RMSE_theta_prop(k_idx) = sqrt(sq_err_theta / (Q_MC * L));
    RMSE_phi_prop(k_idx)   = sqrt(sq_err_phi / (Q_MC * L));
    RMSE_range_prop(k_idx) = sqrt(sq_err_range / (Q_MC * L));
    
    fprintf('Snapshots K = %4d | Theta: %6.4f (CRLB: %6.4f) | Phi: %6.4f (CRLB: %6.4f) | Range: %6.4f (CRLB: %6.4f)\n', ...
            K_current, RMSE_theta_prop(k_idx), CRLB_theta_k(k_idx), RMSE_phi_prop(k_idx), CRLB_phi_k(k_idx), RMSE_range_prop(k_idx), CRLB_range_k(k_idx));
end

save('RMSE_theta_prop.mat', 'RMSE_theta_prop');
save('RMSE_phi_prop.mat',   'RMSE_phi_prop'  );
save('RMSE_range_prop.mat', 'RMSE_range_prop');

save('CRLB_theta_k.mat', 'CRLB_theta_k');
save('CRLB_phi_k.mat',   'CRLB_phi_k'  );
save('CRLB_range_k.mat', 'CRLB_range_k');





