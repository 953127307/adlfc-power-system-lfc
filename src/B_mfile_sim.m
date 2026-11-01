function [t, x_cell, u_cell, P_f, P_f_d, governor_log] = ...
    B_mfile_sim(pam, Time, Step_T, bar_d, mode, w, M_hat_init)
% B_mfile_sim — Pure m-file simulation of 4-area power system (B part).
% Supports three control methods via pam.ctrl_method:
%   1 = Proposed ADLFC (backstepping + adaptive M̂ + low-pass filters)
%   2 = ACE-SMC with theta adaptation (sliding mode control)
%   3 = Conventional ACE integral controller

%% Extract parameters
K_p  = pam.K_p;  mu_p = pam.mu_p;  mu_t = pam.mu_t;  mu_g = pam.mu_g;
R    = pam.R;    xi   = pam.xi;    P_d_0 = pam.P_d_0;
omega = pam.omega;  r_adapt = pam.r;  delta_adapt = pam.delta;
x_d = pam.x_d;  x_0_mat = pam.x_0;
P_d_step = pam.P_d_step_schedule;
n_area = 4;

% Controller method
if isfield(pam, 'ctrl_method'), ctrl_method = pam.ctrl_method;
else, ctrl_method = 1; end

beta_line = xi;
neighbors = {2, [1,3], [2,4], 3};

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


%% Pre-allocate
x_cell = cell(n_area, 1);  u_cell = cell(n_area, 1);
for i = 1:n_area
    x_cell{i} = zeros(Nt, 4);  u_cell{i} = zeros(Nt, 1);
end

%% Initial conditions: start from x_0 (equilibrium for P_d(0))
x = x_0_mat;

% ADLFC init: compute hat_alpha_cur from initial state
if ctrl_method == 1 && strcmp(mode, 'verify')
    if isfield(pam, 'hat_alpha_0') && any(pam.hat_alpha_0(:) ~= 0)
        hat_alpha_cur = pam.hat_alpha_0;
    end
    for i = 1:n_area
        z(i,1) = x(i,1) - x_d(i);
        z(i,2) = x(i,2) - hat_alpha_cur(i,1);
        z(i,3) = x(i,3) - hat_alpha_cur(i,2);
        z(i,4) = x(i,4) - hat_alpha_cur(i,3);
    end
end

%% Pre-compute P_d(t): explicit common mixed-load trajectory when supplied.
if isfield(pam, 'P_d_trajectory')
    P_d_t = pam.P_d_trajectory;
    assert(isequal(size(P_d_t), [Nt, n_area]));
    assert(all(isfinite(P_d_t(:))));
else
    P_d_t = zeros(Nt, n_area);
    if ~isempty(P_d_step)
        for step_idx = 1:size(P_d_step, 1)
            idx = find(t >= P_d_step(step_idx, 1), 1, 'first');
            if ~isempty(idx)
                P_d_t(idx:end, :) = repmat(P_d_step(step_idx, 2:end), Nt - idx + 1, 1);
            end
        end
    end
end
P_d_cur = P_d_t(1,:)';

%% State jump setup (disabled for multi-step simulation)
has_state_jump = false;
t_jump = inf;
x_jump = [];

%% Method-specific state init
theta = zeros(n_area, 1);
xi_prev = zeros(n_area, 1);  xi_prev2 = zeros(n_area, 1);  xi_max = zeros(n_area, 1);
ACE_int = zeros(n_area, 1);
u_int = zeros(n_area, 1);  % integrator state for Case 2 (u = ∫w dt)

%% ADLFC-specific init
M_hat = cell(n_area, 4);
for i = 1:n_area
    M_hat{i,1} = 0;
    for j = 2:4
        if nargin >= 7 && ~isempty(M_hat_init) && ~isempty(M_hat_init{i,j})
            M_hat{i,j} = M_hat_init{i,j};
        else
            M_hat{i,j} = zeros(1, max(1, length(w{i,j}) - 1));
        end
    end
end
z = zeros(n_area, 4);  hat_alpha_cur = zeros(n_area, 3);

if strcmp(mode, 'verify') && ctrl_method == 1
    if isfield(pam, 'hat_alpha_0') && any(pam.hat_alpha_0(:) ~= 0)
        hat_alpha_cur = pam.hat_alpha_0;
    end
    for i = 1:n_area
        z(i,1) = x(i,1) - x_d(i);
        z(i,2) = x(i,2) - hat_alpha_cur(i,1);
        z(i,3) = x(i,3) - hat_alpha_cur(i,2);
        z(i,4) = x(i,4) - hat_alpha_cur(i,3);
    end
end

%% Pre-generate disturbances
rng(1);
d_all = bar_d * (2*rand(Nt, n_area, 4) - 1);

%% ======================== MAIN SIMULATION LOOP ========================
for k = 1:Nt
    % --- State jump injection ---
    if has_state_jump && k > 1 && t(k-1) < t_jump && t(k) >= t_jump
        x = x_jump;
        if ctrl_method == 1
            if isfield(pam, 'hat_alpha_0_post_jump') && any(pam.hat_alpha_0_post_jump(:) ~= 0)
                hat_alpha_cur = pam.hat_alpha_0_post_jump;
            end
            for i = 1:n_area
                z(i,1) = x(i,1) - x_d(i);
                z(i,2) = x(i,2) - hat_alpha_cur(i,1);
                z(i,3) = x(i,3) - hat_alpha_cur(i,2);
                z(i,4) = x(i,4) - hat_alpha_cur(i,3);
            end
        end
    end

    % --- Store state ---
    for i = 1:n_area, x_cell{i}(k,:) = x(i,:); end

    % --- Controller ---
    dot_hat_alpha_cur = zeros(n_area, 3);
    psi_all = cell(n_area, 4);
    u = zeros(n_area, 1);

    if strcmp(mode, 'data')
        tk = t(k);
        for i = 1:n_area
            hat_alpha_cur(i,1) = sin(tk);  hat_alpha_cur(i,2) = sin(tk);
            hat_alpha_cur(i,3) = sin(tk);
            dot_hat_alpha_cur(i,1) = cos(tk);  dot_hat_alpha_cur(i,2) = cos(tk);
            dot_hat_alpha_cur(i,3) = cos(tk);
            u(i) = tk;
            z(i,1) = x(i,1) - x_d(i);
            z(i,2) = x(i,2) - hat_alpha_cur(i,1);
            z(i,3) = x(i,3) - hat_alpha_cur(i,2);
            z(i,4) = x(i,4) - hat_alpha_cur(i,3);
            psi_all{i,1} = z(i,1);
            psi_all{i,2} = build_psi_2(i, x, z, dot_hat_alpha_cur(i,1), neighbors, n_area);
            psi_all{i,3} = [z(i,3); x(i,3); dot_hat_alpha_cur(i,2)];
            psi_all{i,4} = [z(i,4); x(i,2); x(i,4); dot_hat_alpha_cur(i,3)];
        end
    elseif ctrl_method == 1
        % === Method 1: ADLFC (backstepping + adaptive M̂) ===
        for i = 1:n_area
            psi_all{i,1} = z(i,1);
            alpha_1 = w{i,1}' * psi_all{i,1};
            dot_hat_alpha_cur(i,1) = (alpha_1 - hat_alpha_cur(i,1)) / omega(i,1);

            psi_all{i,2} = build_psi_2(i, x, z, dot_hat_alpha_cur(i,1), neighbors, n_area);
            sigma_2 = psi_all{i,2}(2:end);
            alpha_2 = w{i,2}' * psi_all{i,2} - M_hat{i,2} * sigma_2;
            dot_hat_alpha_cur(i,2) = (alpha_2 - hat_alpha_cur(i,2)) / omega(i,2);

            psi_all{i,3} = [z(i,3); x(i,3); dot_hat_alpha_cur(i,2)];
            sigma_3 = psi_all{i,3}(2:end);
            alpha_3 = w{i,3}' * psi_all{i,3} - M_hat{i,3} * sigma_3;
            dot_hat_alpha_cur(i,3) = (alpha_3 - hat_alpha_cur(i,3)) / omega(i,3);

            psi_all{i,4} = [z(i,4); x(i,2); x(i,4); dot_hat_alpha_cur(i,3)];
            sigma_4 = psi_all{i,4}(2:end);
            u_raw = w{i,4}' * psi_all{i,4} - M_hat{i,4} * sigma_4;
            u(i) = u_raw;  % No control-input amplitude saturation.
        end
    elseif ctrl_method == 2
        % === Method 2: ACE-SMC with theta adaptation ===
        % ACE_a = B * (P_f - P_f_d)
        P_f_now = zeros(3,1);
        for l = 1:3
            P_f_now(l) = beta_line(l) * sin(x(l,1) - x(l+1,1));
        end
        ACE_a = pam.B * (P_f_now - [0.015; 0.0125; 0.01]);
        T_theta = 0.1;  W_smc = pam.W;
        for i = 1:n_area
            M1 = 2;  M2 = 0.1;  M3 = 0.01;  M4 = -(M2 + M3);
            M5 = 0.0217 * M1 * K_p(i) / mu_p(i);
            xi_cur = M1*x(i,2) + M2*x(i,3) + M3*x(i,4) + M4*theta(i) + M5*ACE_a(i);
            % SMC law: w = -W * sign(xi - xi_max/2)
            w_smc = -W_smc * sign(xi_cur - 0.5*xi_max(i));
            % Integrate w to produce u (original has Integrator after sign)
            u_int(i) = u_int(i) + w_smc * Step_T;
            u(i) = u_int(i);
            % Peak-hold xi_max (Switch + UnitDelay in original)
            xi_p = xi_prev(i);  xi_p2 = xi_prev2(i);
            if xi_cur <= xi_p && xi_p <= xi_p2
                xi_max(i) = xi_p;  % capture local maximum
            end
            % Theta adaptation: dot_theta = -theta/T + x3/T
            theta(i) = theta(i) + (-theta(i)/T_theta + x(i,3)/T_theta) * Step_T;
            xi_prev2(i) = xi_p;  xi_prev(i) = xi_cur;
        end
    elseif ctrl_method == 3
        % === Method 3: Integral-only ACE (per paper) ===
        P_f_now = zeros(3,1);
        for l = 1:3
            P_f_now(l) = beta_line(l) * sin(x(l,1) - x(l+1,1));
        end
        K_I3 = 0.2;  B_mat = pam.B;  P_f_d_ref = [0.015; 0.0125; 0.01];
        for i = 1:n_area
            b_i = 1/R(i) + 1/K_p(i);
            ACE_i = b_i * x(i,2);
            for l = 1:3
                ACE_i = ACE_i + B_mat(i,l) * (P_f_now(l) - P_f_d_ref(l));
            end
            ACE_int(i) = ACE_int(i) + ACE_i * Step_T;
            u_raw = -K_I3 * ACE_int(i);
            u(i) = u_raw;
        end
    end

    for i = 1:n_area, u_cell{i}(k) = u(i); end

    % --- Plant dynamics ---
    dx = zeros(n_area, 4);
    for i = 1:n_area
        coupling = 0;
        for nb = neighbors{i}
            line_idx = min(i, nb);
            coupling = coupling - beta_line(line_idx) * sin(x(i,1) - x(nb,1));
        end
        dx(i,1) = x(i,2) + d_all(k,i,1);
        dx(i,2) = -(1/mu_p(i))*x(i,2) + (K_p(i)/mu_p(i))*(x(i,3) - P_d_cur(i) + coupling) + d_all(k,i,2);
        dx(i,3) = (1/mu_t(i))*x(i,4) - (1/mu_t(i))*x(i,3) + d_all(k,i,3);
        dx(i,4) = (1/mu_g(i))*u(i) - (1/(mu_g(i)*R(i)))*x(i,2) - (1/mu_g(i))*x(i,4) + d_all(k,i,4);
        governor_log.raw_rate(k,i) = dx(i,4);
        dx(i,4) = max(-governor_rate_limit, min(governor_rate_limit, dx(i,4)));
        governor_log.actual_rate(k,i) = dx(i,4);
    end

    % --- Update ---
    if k < Nt
        if strcmp(mode, 'verify') && ctrl_method == 1
            for i = 1:n_area
                for j = 2:4
                    sigma_ij = psi_all{i,j}(2:end);
                    M_hat{i,j} = M_hat{i,j} + r_adapt(i,j) * (z(i,j)*sigma_ij' - delta_adapt(i,j)*M_hat{i,j}) * Step_T;
                end
            end
        end
        x = x + dx * Step_T;
        P_d_cur = P_d_t(min(k+1, Nt), :)';
        if ctrl_method == 1
            hat_alpha_cur = hat_alpha_cur + dot_hat_alpha_cur * Step_T;
            for i = 1:n_area
                z(i,1) = x(i,1) - x_d(i);
                z(i,2) = x(i,2) - hat_alpha_cur(i,1);
                z(i,3) = x(i,3) - hat_alpha_cur(i,2);
                z(i,4) = x(i,4) - hat_alpha_cur(i,3);
            end
        end
    end
end

%% --- Tie-line power ---
P_f = zeros(Nt, 3);
for k = 1:Nt
    P_f(k,1) = beta_line(1) * sin(x_cell{1}(k,1) - x_cell{2}(k,1));
    P_f(k,2) = beta_line(2) * sin(x_cell{2}(k,1) - x_cell{3}(k,1));
    P_f(k,3) = beta_line(3) * sin(x_cell{3}(k,1) - x_cell{4}(k,1));
end
P_f_d = repmat([0.015, 0.0125, 0.01], Nt, 1);
end

%% ======================== HELPER ============================
function psi_2 = build_psi_2(i, x, z, dot_hat_alpha_1, neighbors, n_area)
nb = neighbors{i};
sin_terms = zeros(length(nb), 1);
for n = 1:length(nb)
    sin_terms(n) = sin(x(i,1) - x(nb(n),1));
end
psi_2 = [z(i,2); sin_terms; x(i,2); dot_hat_alpha_1; 1];
end
