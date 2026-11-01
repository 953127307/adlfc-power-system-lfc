function [t, x_cell, u_cell, P_f, P_f_d, dot_z_cell, psi_cell, d_cell, ...
          hat_alpha, hat_M_out, governor_log] = ...
    A_mfile_sim(pam, Time, Step_T, bar_d, mode, w, M_hat_init)
% A_mfile_sim — Pure m-file simulation of 4-area ADLFC power system.
% Replaces both A_S1.slx (data collection) and A_S3.slx (verification).
%
% Inputs:
%   pam         — parameter struct (from P_sim.m)
%   Time        — simulation duration (s)
%   Step_T      — integration step size (s)
%   bar_d       — disturbance bound (scalar)
%   mode        — 'data'  = data collection with sin(t) controller
%                 'verify'= closed-loop with backstepping controller
%   w           — 4x4 cell of controller weight vectors (only for 'verify')
%   M_hat_init  — 4x4 cell of initial adaptive params (only for 'verify')
%
% Outputs:
%   t           — time vector (Nt x 1)
%   x_cell      — 4x1 cell, x_cell{i} = state matrix (Nt x 4)
%   u_cell      — 4x1 cell, u_cell{i} = control input vector (Nt x 1)
%   P_f         — tie-line power (Nt x 3): [P12, P23, P34]
%   P_f_d       — desired tie-line power (Nt x 3)
%   dot_z_cell  — 4x4 cell, dot_z_cell{i,j} = error derivatives (Nt x 1)
%   psi_cell    — 4x4 cell, psi_cell{i,j} = basis functions (Nt x m_ij)
%   d_cell      — 4x4 cell, d_cell{i,j} = disturbances (Nt x 1)
%   hat_alpha   — Nt x 4 x 3 array of filtered virtual controls
%   hat_M_out   — 4x4 cell of adaptive parameter trajectories

%% Extract parameters
K_p  = pam.K_p;
mu_p = pam.mu_p;
mu_t = pam.mu_t;
mu_g = pam.mu_g;
R    = pam.R;
xi   = pam.xi;
P_d_0 = pam.P_d_0;
omega = pam.omega;   % 4x3 filter time constants
r_adapt = pam.r;     % 4x4 adaptive gains
delta_adapt = pam.delta;  % 4x4 sigma-modification gains

x_d = pam.x_d;       % 4x1 reference angles
x_0_mat = pam.x_0;   % 4x4 initial states (each row = one area)

P_d_step = pam.P_d_step_schedule;  % load step schedule
if isfield(pam, 'Pg_step_schedule')
    Pg_step = pam.Pg_step_schedule;    % generator state jump schedule
else
    Pg_step = [];
end

n_area = 4;

%% Derived parameters
beta_line = xi;  % [β12, β23, β34] — line synchronizing coefficients

%% Neighbor sets for each area (chain topology: 1-2-3-4)
neighbors = cell(n_area, 1);
neighbors{1} = 2;
neighbors{2} = [1, 3];
neighbors{3} = [2, 4];
neighbors{4} = 3;

%% Determine simulation steps
Nt = round(Time / Step_T) + 1;
t = (0:Nt-1)' * Step_T;

% R2 Comment 2.2: rate limit only; no added governor/turbine amplitude limit.
% Clamp the complete governor derivative, including the injected d_4 term.
governor_rate_limit = inf;
if isfield(pam, 'governor_rate_limit')
    governor_rate_limit = pam.governor_rate_limit;
end
assert(isscalar(governor_rate_limit) && governor_rate_limit > 0);
governor_log.raw_rate = zeros(Nt, n_area);
governor_log.actual_rate = zeros(Nt, n_area);
governor_log.limit = governor_rate_limit;


%% Pre-allocate storage
x_cell = cell(n_area, 1);
u_cell = cell(n_area, 1);
for i = 1:n_area
    x_cell{i} = zeros(Nt, 4);
    u_cell{i} = zeros(Nt, 1);
end

dot_z_cell = cell(n_area, 4);
psi_cell = cell(n_area, 4);
d_cell = cell(n_area, 4);
for i = 1:n_area
    for j = 1:4
        dot_z_cell{i,j} = zeros(Nt, 1);
    end
end

hat_alpha = zeros(Nt, n_area, 3);
hat_M_out = cell(n_area, 4);

%% Initialize M_hat (adaptive parameters) for verification mode
M_hat = cell(n_area, 4);
if strcmp(mode, 'verify')
    for i = 1:n_area
        M_hat{i,1} = 0;
        hat_M_out{i,1} = zeros(Nt, 1);
        for j = 2:4
            m_ij = length(w{i,j});
            if nargin >= 7 && ~isempty(M_hat_init) && ~isempty(M_hat_init{i,j})
                M_hat{i,j} = M_hat_init{i,j};
            else
                M_hat{i,j} = zeros(1, m_ij - 1);
            end
            hat_M_out{i,j} = zeros(Nt, m_ij - 1);
        end
    end
else
    % Data collection mode: no adaptive terms needed
    for i = 1:n_area
        for j = 1:4
            M_hat{i,j} = 0;
            hat_M_out{i,j} = zeros(Nt, 1);
        end
    end
end

%% Initial conditions
x = x_0_mat;  % n_area x 4

%% Pre-compute P_d(t) for load steps
P_d_t = zeros(Nt, n_area);
if ~isempty(P_d_step)
    for step_idx = 1:size(P_d_step, 1)
        t_step = P_d_step(step_idx, 1);
        idx = find(t >= t_step, 1, 'first');
        if ~isempty(idx)
            P_d_t(idx:end, :) = repmat(P_d_step(step_idx, 2:end), Nt - idx + 1, 1);
        end
    end
else
    P_d_t = repmat(P_d_0', Nt, 1);
end
P_d_cur = P_d_t(1,:)';

%% Initialize errors and filter states
z = zeros(n_area, 4);
hat_alpha_cur = zeros(n_area, 3);

% Compute initial z and hat_alpha for both modes
for i = 1:n_area
    z(i,1) = x(i,1) - x_d(i);
end

if strcmp(mode, 'verify')
    % Use hat_alpha_0 if available, otherwise compute from weights
    if isfield(pam, 'hat_alpha_0') && ~isempty(pam.hat_alpha_0) ...
            && any(pam.hat_alpha_0(:) ~= 0)
        hat_alpha_cur = pam.hat_alpha_0;
    else
        for i = 1:n_area
            alpha_1 = w{i,1}' * z(i,1);
            hat_alpha_cur(i,1) = alpha_1;
            z(i,2) = x(i,2) - hat_alpha_cur(i,1);
            dot_ha_1 = (alpha_1 - hat_alpha_cur(i,1)) / omega(i,1);

            psi_2 = build_psi_2(i, x, z, dot_ha_1, neighbors);
            sigma_2 = psi_2(2:end);
            alpha_2 = w{i,2}' * psi_2 - M_hat{i,2} * sigma_2;
            hat_alpha_cur(i,2) = alpha_2;
            z(i,3) = x(i,3) - hat_alpha_cur(i,2);
            dot_ha_2 = (alpha_2 - hat_alpha_cur(i,2)) / omega(i,2);

            psi_3 = [z(i,3); x(i,3); dot_ha_2];
            sigma_3 = psi_3(2:end);
            alpha_3 = w{i,3}' * psi_3 - M_hat{i,3} * sigma_3;
            hat_alpha_cur(i,3) = alpha_3;
        end
    end
    % Compute remaining z
    for i = 1:n_area
        z(i,1) = x(i,1) - x_d(i);
        z(i,2) = x(i,2) - hat_alpha_cur(i,1);
        z(i,3) = x(i,3) - hat_alpha_cur(i,2);
        z(i,4) = x(i,4) - hat_alpha_cur(i,3);
    end
else
    % Data mode: α̂(0)=sin(0)=0, u(0)=0
    for i = 1:n_area
        z(i,1) = x(i,1) - x_d(i);
        z(i,2) = x(i,2) - 0;
        z(i,3) = x(i,3) - 0;
        z(i,4) = x(i,4) - 0;
    end
end

%% Seed random generator for reproducibility
rng(1);

% Pre-generate all disturbances: d_{i,j}(t) ~ Uniform[-bar_d, bar_d]
d_all = zeros(Nt, n_area, 4);
for i = 1:n_area
    for j = 1:4
        d_all(:,i,j) = bar_d * (2*rand(Nt,1) - 1);
    end
end

%% ======================== MAIN SIMULATION LOOP ========================
for k = 1:Nt
    % --- Store current state ---
    for i = 1:n_area
        x_cell{i}(k,:) = x(i,:);
    end

    % --- Controller computation ---
    dot_hat_alpha_cur = zeros(n_area, 3);
    psi_all = cell(n_area, 4);
    u = zeros(n_area, 1);

    if strcmp(mode, 'data')
        % Data collection per manuscript Section 4.1:
        %   α̂_{i,l}(t) = sin(t), dot_α̂_{i,l}(t) = cos(t)
        %   u_i(t) = t
        % These are DIRECTLY SET (filter dynamics bypassed).
        tk = t(k);
        for i = 1:n_area
            hat_alpha_cur(i,1) = sin(tk);
            hat_alpha_cur(i,2) = sin(tk);
            hat_alpha_cur(i,3) = sin(tk);
            dot_hat_alpha_cur(i,1) = cos(tk);
            dot_hat_alpha_cur(i,2) = cos(tk);
            dot_hat_alpha_cur(i,3) = cos(tk);

            u(i) = tk;  % u_j = t

            % Compute z and psi from known hat_alpha
            z(i,1) = x(i,1) - x_d(i);
            z(i,2) = x(i,2) - hat_alpha_cur(i,1);
            z(i,3) = x(i,3) - hat_alpha_cur(i,2);
            z(i,4) = x(i,4) - hat_alpha_cur(i,3);

            psi_all{i,1} = z(i,1);
            psi_all{i,2} = build_psi_2(i, x, z, dot_hat_alpha_cur(i,1), neighbors);
            psi_all{i,3} = [z(i,3); x(i,3); dot_hat_alpha_cur(i,2)];
            psi_all{i,4} = [z(i,4); x(i,2); x(i,4); dot_hat_alpha_cur(i,3)];
        end
    else
        % Verification: full backstepping controller with adaptive M̂
        for i = 1:n_area
            % ---- Step 1 ----
            psi_all{i,1} = z(i,1);
            alpha_1 = w{i,1}' * psi_all{i,1};
            dot_hat_alpha_cur(i,1) = (alpha_1 - hat_alpha_cur(i,1)) / omega(i,1);

            % ---- Step 2 ----
            psi_all{i,2} = build_psi_2(i, x, z, dot_hat_alpha_cur(i,1), neighbors);
            sigma_2 = psi_all{i,2}(2:end);
            alpha_2 = w{i,2}' * psi_all{i,2} - M_hat{i,2} * sigma_2;
            dot_hat_alpha_cur(i,2) = (alpha_2 - hat_alpha_cur(i,2)) / omega(i,2);

            % ---- Step 3 ----
            psi_all{i,3} = [z(i,3); x(i,3); dot_hat_alpha_cur(i,2)];
            sigma_3 = psi_all{i,3}(2:end);
            alpha_3 = w{i,3}' * psi_all{i,3} - M_hat{i,3} * sigma_3;
            dot_hat_alpha_cur(i,3) = (alpha_3 - hat_alpha_cur(i,3)) / omega(i,3);

            % ---- Step 4 (actual control) ----
            psi_all{i,4} = [z(i,4); x(i,2); x(i,4); dot_hat_alpha_cur(i,3)];
            sigma_4 = psi_all{i,4}(2:end);
            u_raw = w{i,4}' * psi_all{i,4} - M_hat{i,4} * sigma_4;
            u(i) = u_raw;  % No control-input amplitude saturation.
        end
    end

    % Store control
    for i = 1:n_area
        u_cell{i}(k) = u(i);
    end

    % --- Plant dynamics ---
    dx = zeros(n_area, 4);
    for i = 1:n_area
        % Coupling term: -Σ β_{ik} * sin(x_i - x_k)
        coupling = 0;
        for nb = neighbors{i}
            line_idx = min(i, nb);  % maps (1,2)→1, (2,3)→2, (3,4)→3
            coupling = coupling - beta_line(line_idx) * sin(x(i,1) - x(nb,1));
        end

        % State equations from manuscript Eq. (2)
        dx(i,1) = x(i,2) + d_all(k,i,1);
        dx(i,2) = -(1/mu_p(i))*x(i,2) ...
                + (K_p(i)/mu_p(i)) * (x(i,3) - P_d_cur(i) + coupling) ...
                + d_all(k,i,2);
        dx(i,3) = (1/mu_t(i))*x(i,4) - (1/mu_t(i))*x(i,3) + d_all(k,i,3);
        dx(i,4) = (1/mu_g(i))*u(i) ...
                - (1/(mu_g(i)*R(i)))*x(i,2) - (1/mu_g(i))*x(i,4) ...
                + d_all(k,i,4);
        governor_log.raw_rate(k,i) = dx(i,4);
        dx(i,4) = max(-governor_rate_limit, min(governor_rate_limit, dx(i,4)));
        governor_log.actual_rate(k,i) = dx(i,4);
    end

    % --- Compute dot_z (error derivatives) ---
    for i = 1:n_area
        dot_z_cell{i,1}(k) = dx(i,1);  % dot_z_{i,1} = dot_x_{i,1}
        dot_z_cell{i,2}(k) = dx(i,2) - dot_hat_alpha_cur(i,1);
        dot_z_cell{i,3}(k) = dx(i,3) - dot_hat_alpha_cur(i,2);
        dot_z_cell{i,4}(k) = dx(i,4) - dot_hat_alpha_cur(i,3);
    end

    % --- Store basis functions ---
    for i = 1:n_area
        for j = 1:4
            if isempty(psi_cell{i,j})
                psi_cell{i,j} = zeros(Nt, length(psi_all{i,j}));
            end
            psi_cell{i,j}(k,:) = psi_all{i,j}';
        end
    end

    % --- Store disturbances ---
    for i = 1:n_area
        for j = 1:4
            if isempty(d_cell{i,j})
                d_cell{i,j} = zeros(Nt, 1);
            end
            d_cell{i,j}(k) = d_all(k,i,j);
        end
    end

    % --- Store hat_alpha and hat_M ---
    hat_alpha(k,:,:) = hat_alpha_cur;
    for i = 1:n_area
        for j = 2:4
            if isempty(hat_M_out{i,j})
                hat_M_out{i,j} = zeros(Nt, length(M_hat{i,j}));
            end
            hat_M_out{i,j}(k,:) = M_hat{i,j};
        end
        if isempty(hat_M_out{i,1})
            hat_M_out{i,1} = zeros(Nt, 1);
        end
    end

    % --- Update for next step ---
    if k < Nt
        % Update M̂ (adaptive law from manuscript Eq. 3:hat_M)
        if strcmp(mode, 'verify')
            for i = 1:n_area
                for j = 2:4
                    sigma_ij = psi_all{i,j}(2:end);
                    dot_M = r_adapt(i,j) * (z(i,j) * sigma_ij' - delta_adapt(i,j) * M_hat{i,j});
                    M_hat{i,j} = M_hat{i,j} + dot_M * Step_T;
                end
            end
        end

        % Forward Euler integration
        x = x + dx * Step_T;

        % Apply generator state jumps at scheduled times
        if ~isempty(Pg_step)
            for step_idx = 1:size(Pg_step, 1)
                t_pg = Pg_step(step_idx, 1);
                if abs(t(k+1) - t_pg) < Step_T/2
                    for i = 1:n_area
                        if Pg_step(step_idx, i+1) ~= 0
                            x(i,3) = x(i,3) + Pg_step(step_idx, i+1);
                        end
                    end
                end
            end
        end

        % Update P_d
        P_d_cur = P_d_t(min(k+1, Nt), :)';

        % Update filter states (hat_alpha)
        hat_alpha_cur = hat_alpha_cur + dot_hat_alpha_cur * Step_T;

        % Update errors (z)
        for i = 1:n_area
            z(i,1) = x(i,1) - x_d(i);
            z(i,2) = x(i,2) - hat_alpha_cur(i,1);
            z(i,3) = x(i,3) - hat_alpha_cur(i,2);
            z(i,4) = x(i,4) - hat_alpha_cur(i,3);
        end
    end
end

%% --- Compute tie-line power ---
P_f = zeros(Nt, 3);
for k = 1:Nt
    P_f(k,1) = beta_line(1) * sin(x_cell{1}(k,1) - x_cell{2}(k,1));  % P12
    P_f(k,2) = beta_line(2) * sin(x_cell{2}(k,1) - x_cell{3}(k,1));  % P23
    P_f(k,3) = beta_line(3) * sin(x_cell{3}(k,1) - x_cell{4}(k,1));  % P34
end

% Desired tie-line power (constant)
P_f_d = repmat([0.015, 0.0125, 0.01], Nt, 1);

end

%% ======================== HELPER FUNCTIONS ============================

function psi_2 = build_psi_2(i, x, z, dot_hat_alpha_1, neighbors)
% Build ψ_{i,2} = [z_{i,2}; sin(x_{i,1}-x_{k,1}) for k∈A_i; x_{i,2}; dot_α̂_{i,1}; 1]

nb = neighbors{i};
n_nb = length(nb);
sin_terms = zeros(n_nb, 1);
for n = 1:n_nb
    sin_terms(n) = sin(x(i,1) - x(nb(n),1));
end

psi_2 = [z(i,2); sin_terms; x(i,2); dot_hat_alpha_1; 1];
end
