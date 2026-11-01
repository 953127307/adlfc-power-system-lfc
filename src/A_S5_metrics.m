function metrics = A_S5_metrics(t, x_cell, P_f, P_f_d, pam, Step_T)
% A_S5_metrics  计算 LFC 全部性能指标
%
% 输入:
%   t       - 时间向量 (Nt x 1)
%   x_cell  - 4x1 cell, x_cell{i} 为区域 i 的状态矩阵 (Nt x 4)
%             列: [rotor_angle, freq_dev, turbine_power, governor_power]
%   P_f     - 联络线功率 (Nt x 3), 列: [P_12, P_23, P_34]
%   P_f_d   - 期望联络线功率 (Nt x 3)
%   pam     - 参数结构体 (含 B, B_f, P_d_step_schedule, K_p, mu_p)
%   Step_T  - 仿真步长
%
% 输出:
%   metrics - 结构体, 含 ACE, RoCoF, Nadir, IAE, ISE, RMS, SettlingTime

Nt = length(t);
n_area = 4;

% 确保为列向量
t = t(:);

% 处理 P_f_d 维度 (可能为 N×3 或 1×3 或 3×1)
if size(P_f_d, 1) == 3 && size(P_f_d, 2) == 1
    P_f_d_ref = P_f_d';  % 3x1 -> 1x3
elseif size(P_f_d, 2) == 3
    P_f_d_ref = P_f_d(1,:);  % 取第一行
else
    P_f_d_ref = [0.015, 0.0125, 0.01];  % fallback
end
% 扩展为与 P_f 相同行数
P_f_d_mat = repmat(P_f_d_ref, size(P_f, 1), 1);

%% 1. 计算 ACE (Area Control Error)
% ACE_i = Delta_P_tie,i + B_f(i) * Delta_f_i
% Delta_P_tie = B * (P_f - P_f_d)
B = pam.B;          % 4x3 关联矩阵
B_f = pam.B_f(:);   % 4x1 频率偏差系数

ACE_cell = cell(n_area, 1);
for i = 1:n_area
    Delta_P_tie_i = (P_f - P_f_d_mat) * B(i,:)';  % Nt x 1
    Delta_f_i = x_cell{i}(:, 2);                   % 频率偏差 gamma_i
    ACE_cell{i} = Delta_P_tie_i + B_f(i) * Delta_f_i;
end

%% 2. 积分指标: IAE, ISE, RMS of ACE
IAE_ACE = zeros(n_area, 1);
ISE_ACE = zeros(n_area, 1);
RMS_ACE = zeros(n_area, 1);

for i = 1:n_area
    ace = ACE_cell{i};
    IAE_ACE(i) = trapz(t, abs(ace));
    ISE_ACE(i) = trapz(t, ace.^2);
    RMS_ACE(i) = sqrt(mean(ace.^2));
end

%% 3. RoCoF 和 Nadir (对每个负荷阶跃事件+发电机事件)
P_d_schedule = pam.P_d_step_schedule;
step_times = P_d_schedule(2:end, 1);  % 负荷阶跃时刻 (排除 t=0 baseline)
% 合并发电机损失事件
if isfield(pam, 'Pg_step_schedule') && ~isempty(pam.Pg_step_schedule)
    pg_times = pam.Pg_step_schedule(:, 1);
    step_times = unique([step_times; pg_times]);
end
step_times = sort(step_times);
n_steps = length(step_times);

RoCoF = zeros(n_area, n_steps);
Nadir = zeros(n_area, n_steps);

% 阶跃后分析窗
window_RoCoF = 0.5;    % RoCoF: 阶跃后 0.5 秒内
window_Nadir = 10;     % Nadir: 阶跃后 10 秒内

for s = 1:n_steps
    t_step = step_times(s);
    idx_start = find(t >= t_step, 1, 'first');
    if isempty(idx_start); continue; end

    % RoCoF 窗
    idx_end_RoCoF = find(t >= t_step + window_RoCoF, 1, 'first');
    if isempty(idx_end_RoCoF); idx_end_RoCoF = Nt; end

    % Nadir 窗
    idx_end_Nadir = find(t >= t_step + window_Nadir, 1, 'first');
    if isempty(idx_end_Nadir); idx_end_Nadir = Nt; end

    for i = 1:n_area
        gamma_i = x_cell{i}(:, 2);
        % RoCoF: 频率变化率的最大绝对值
        d_gamma = diff(gamma_i(idx_start:idx_end_RoCoF)) / Step_T;
        if ~isempty(d_gamma)
            RoCoF(i, s) = max(abs(d_gamma));
        end
        % Nadir: 阶跃后频率最低点
        Nadir(i, s) = min(gamma_i(idx_start:idx_end_Nadir));
    end
end

%% 4. Settling Time (每次阶跃后 ACE 进入稳态带，需持续 1.0s)
% 基准：1% × 当前 Case 所有区域中的最高峰值偏差，四区域共用同一 band
% 参考值：以 0 为稳态目标（LFC 目标即为 ACE→0）
% No 2 s exit-time or twice-band-amplitude screening.
tolerance_frac = 0.01;
hold_t = 1.0;
StepSettling = NaN(n_area, n_steps);
settling_time = zeros(n_area, 1);

win_ends = [step_times(2:end); t(end)];
hold_N = ceil(hold_t / Step_T - 1e-9) + 1; % 1001 samples span 1 s at dt=0.001 s

% 各区域全程峰值偏差（以 0 为参考）
idx_band_start = find(t >= 0, 1, 'first');
if isempty(idx_band_start); idx_band_start = 1; end

for i = 1:n_area
    ace_full = ACE_cell{i};
    peak_dev(i) = max(abs(ace_full(idx_band_start:end)));
    if peak_dev(i) < 1e-12; peak_dev(i) = 1e-12; end
end

% 当前 Case 统一 band = 1% × 全局最高峰值偏差
band = tolerance_frac * max(peak_dev);
band_ref = band * ones(n_area, 1);

for s = 1:n_steps
    t0 = step_times(s);
    t_end = win_ends(s);
    if s < n_steps
        idx_win = find(t >= t0 & t < t_end);
    else
        idx_win = find(t >= t0 & t <= t_end);
    end

    for i = 1:n_area
        ace_win = ACE_cell{i}(idx_win);
        t_win = t(idx_win);
        N_win = length(ace_win);
        if N_win < hold_N
            StepSettling(i,s) = inf;
            continue;
        end
        in_band = abs(ace_win) <= band;
        idx_first_out = find(~in_band, 1, 'first');
        if isempty(idx_first_out)
            StepSettling(i,s) = 0;
            continue;
        end
        % No latency or excursion-amplitude screening is applied.
        % Find the first re-entry followed by at least one full second in band.
        settled_idx = 0;
        cnt = 0;
        for k = idx_first_out:N_win
            if in_band(k)
                cnt = cnt + 1;
                if cnt >= hold_N
                    settled_idx = k - hold_N + 1;
                    break;
                end
            else
                cnt = 0;
            end
        end
        if settled_idx > 0
            StepSettling(i,s) = t_win(settled_idx) - t0;
        else
            StepSettling(i,s) = inf;
        end
    end
end
settling_time = max(StepSettling, [], 2);
settling_time(isnan(settling_time)) = 0;

%% 5. 汇总输出
metrics.ACE = ACE_cell;
metrics.IAE_ACE = IAE_ACE;
metrics.ISE_ACE = ISE_ACE;
metrics.RMS_ACE = RMS_ACE;
metrics.RoCoF = RoCoF;
metrics.Nadir = Nadir;
metrics.StepSettling = StepSettling;
metrics.SettlingTime = settling_time;
metrics.step_times = step_times;
metrics.band = band;
metrics.band_ref = band_ref;

%% 6. 打印摘要
fprintf('\n========== Band Reference ==========\n');
fprintf('%-12s %12s\n', 'Area', 'peak_dev');
fprintf('------------------------------------\n');
for i = 1:n_area
    fprintf('Area %d     %12.6e\n', i, peak_dev(i));
end
fprintf('Common band = %.6e (1%% x max peak_dev, reference = 0)\n\n', band);
fprintf('========== LFC Performance Metrics ==========\n');
fprintf('%-12s %10s %10s %10s %10s\n', 'Area', 'IAE_ACE', 'ISE_ACE', 'RMS_ACE', 'SetlTime');
fprintf('----------------------------------------------\n');
for i = 1:n_area
    fprintf('Area %d     %10.4f %10.6f %10.6f %10.2f\n', ...
        i, IAE_ACE(i), ISE_ACE(i), RMS_ACE(i), settling_time(i));
end
fprintf('\nRoCoF (rad/s^2) at each load step:\n');
fprintf('%-12s', 'Area');
for s = 1:n_steps
    fprintf(' t=%-5.0fs', step_times(s));
end
fprintf('\n');
for i = 1:n_area
    fprintf('Area %d     ', i);
    for s = 1:n_steps
        fprintf('%10.6f', RoCoF(i,s));
    end
    fprintf('\n');
end
fprintf('\nNadir (rad/s) at each load step:\n');
fprintf('%-12s', 'Area');
for s = 1:n_steps
    fprintf(' t=%-5.0fs', step_times(s));
end
fprintf('\n');
for i = 1:n_area
    fprintf('Area %d     ', i);
    for s = 1:n_steps
        fprintf('%10.6f', Nadir(i,s));
    end
    fprintf('\n');
end
fprintf('\nSettling Time (s) at each load step:\n');
fprintf('%-12s', 'Area');
for s = 1:n_steps
    fprintf(' t=%-5.0fs', step_times(s));
end
fprintf('\n');
for i = 1:n_area
    fprintf('Area %d     ', i);
    for s = 1:n_steps
        if isnan(StepSettling(i,s))
            fprintf('%10s', 'N/A');
        elseif isinf(StepSettling(i,s))
            fprintf('%10s', 'Inf');
        else
            fprintf('%10.2f', StepSettling(i,s));
        end
    end
    fprintf('\n');
end
fprintf('==============================================\n');

end
