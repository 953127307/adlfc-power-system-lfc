% A_S1_data_collection.m
% Data collection per manuscript Section 4.1:
%   α̂_{i,l}(t) = sin(t),  u_i(t) = t
%   x_{i,j}(0) = 0,  P_{d_i}(0) applied
%   Simulation: 1s, sampling interval T = 0.1s

fprintf('=== Data Collection (m-file: alpha_hat=sin(t), u=t, %.0fs) ===\n', Time_data);

% Run m-file simulation in data collection mode at Step_T resolution
[t_data, x_cell, u_cell, ~, ~, dot_z_cell, psi_cell, d_cell] = ...
    A_mfile_sim(pam, Time_data, Step_T, bar_d, 'data', [], []);

% Decimate to sampling interval T = 0.1s (per manuscript: "sampling interval T=0.1")
dec_factor = round(T / Step_T);  % 0.1 / 0.001 = 100
idx_sample = 1:dec_factor:length(t_data);
fprintf('  Collected %d steps, decimated to %d samples (T=%.1fs).\n', ...
    length(t_data), length(idx_sample), T);

%% Build data matrices (decimated to T=0.1 sampling per manuscript)
% State data
X_i = cell(n_1, 1);
for i = 1:n_1
    X_i{i} = x_cell{i}(idx_sample, :)';  % 4 x K
end

K_dec = length(idx_sample);
X = cell(n_1, n_2);
for i = 1:4
    for j = 1:4
        X{i,j} = X_i{i}(j, :);  % 1 x K
    end
end

% Control data
U = cell(n_1, 1);
for i = 1:n_1
    U{i} = u_cell{i}(idx_sample)';  % 1 x K
end

% Error derivative data
Z_i = cell(n_1, 1);
for i = 1:n_1
    Z_i{i} = zeros(4, K_dec);
    for j = 1:4
        Z_i{i}(j, :) = dot_z_cell{i,j}(idx_sample)';  % 1 x K
    end
end

Z = cell(n_1, n_2);
for i = 1:4
    for j = 1:4
        Z{i,j} = Z_i{i}(j, :);  % 1 x K
    end
end

% Basis function data
Psi = cell(n_1, n_2);
for i = 1:n_1
    for j = 1:n_2
        Psi{i,j} = psi_cell{i,j}(idx_sample, :)';  % m_{i,j} x K
    end
end

% Disturbance data (for verification only, NOT used in controller design)
D_i = cell(n_1, 1);
for i = 1:n_1
    D_i{i} = zeros(4, K_dec);
    for j = 1:4
        D_i{i}(j, :) = d_cell{i,j}(idx_sample)';  % 1 x K
    end
end

D = cell(n_1, n_2);
for i = 1:4
    for j = 1:4
        D{i,j} = D_i{i}(j, :);  % 1 x K
    end
end

fprintf('  Data matrices ready. K = %d samples.\n', size(U{1}, 2));
