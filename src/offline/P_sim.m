rng(1)

% 系统维度
n_1 = 4;
n_2 = 4;

% 数据采集
Time_data = 1;
T = 0.1;

% 实际仿真
Time = 600;

% 仿真步长
Step_T = 0.001;

% 状态初始值
x_d = zeros(n_1,1);
x_d(1) = 0.4;
x_d(2) = x_d(1)-asin(0.15);
x_d(3) = x_d(2)-asin(0.125);
x_d(4) = x_d(3)-asin(0.1);

pam.x_d = x_d;

x_0 = [0,0,0,0;
       0,0,0,0;
       0,0,0,0;
       0,0,0,0];

pam.x_0 = x_0;

% 系统参数
pam.K_p=[120;112.5;115;118.5];
pam.mu_p=[21;25;23;22];
pam.P_d_0=[0.01;0.014;0.012;0.013];
pam.mu_t=[0.3;0.33;0.35;0.28];
pam.mu_g=[0.08;0.072;0.07;0.081];
pam.R=[2.5;2.7;2.6;2.8];
pam.xi=[0.1;0.1;0.1];

eta = [ones(4,1),(pam.K_p./pam.mu_p),(1./pam.mu_t),(1./pam.mu_g)];

pam.rho_11 = 0;
pam.rho_12 = [-pam.xi(1)*pam.K_p(1)/pam.mu_p(1);-1/pam.mu_p(1);-1;-pam.P_d_0(1)*pam.K_p(1)/pam.mu_p(1)];
pam.rho_13 = [-1/pam.mu_t(1);-1;0];
pam.rho_14 = [-1/(pam.mu_g(1)*pam.R(1));-1/pam.mu_g(1);-1;0];
pam.rho_21 = 0;
pam.rho_22 = [pam.xi(1)*pam.K_p(2)/pam.mu_p(2);pam.xi(2)*pam.K_p(2)/pam.mu_p(2);-1/pam.mu_p(2);-1;-pam.P_d_0(2)*pam.K_p(2)/pam.mu_p(2)];
pam.rho_23 = [-1/pam.mu_t(2);-1;0];
pam.rho_24 = [-1/(pam.mu_g(2)*pam.R(2));-1/pam.mu_g(2);-1;0];

bar_d_case = [1e-10; 0.05; 0.1];  % Case1: 无噪声; Case2: 中等; Case3: 较强

% 负荷阶跃调度: [时间(s), P_d_Area1, P_d_Area2, P_d_Area3, P_d_Area4]
P_d_step_schedule = [
      0,  0.020, 0.028, 0.024, 0.026;   % t=0s baseline
     60,  0.0325,0.028, 0.024, 0.026;   % Area1 负荷 +0.0125
    200,  0.0325,0.0405,0.024, 0.026;   % Area2 负荷 +0.0125
    360,  0.0325,0.0405,0.0365,0.026;   % Area3 负荷 +0.0125
    500,  0.020, 0.028, 0.0365,0.026;   % 恢复至baseline，Area3不恢复
];

pam.B_f = 1 ./ pam.R + 1 ./ pam.K_p;   % 频率偏差系数，用于ACE计算 (标准定义)
pam.P_d_step_schedule = P_d_step_schedule;

% 控制器参数
karpi = [2,520,3,-120;
         2,840,3,-150;
         3,1100,3,-150;
         16,520,3,-120];

pam.omega = 0.2*[0.05,0.05,0.1;
             0.05,0.05,0.1;
             0.05,0.05,0.1;
             0.05,0.05,0.1];

pam.r = 1*[0,0.5,0.5,10;
           0,0.5,0.5,10;
           0,0.5,0.5,10;
           0,0.5,0.5,10];

pam.delta = 1*[0,0.01,0.001,0.01;
               0,0.01,0.001,0.01;
               0,0.01,0.001,0.01;
               0,0.01,0.001,0.01];

pam.B = [ 1  0  0;
     -1  1  0;
      0 -1  1;
      0  0 -1];

pam.Xi = diag([pam.xi(1),pam.xi(2),pam.xi(3)]);

P_f_d = [0.015;0.0125;0.01];
% Preserve the original virtual feedback gains; choose kappa from data.
pam.feedback_gain = ones(4,4);

% R2 Comment 2.2: governor output rate, p.u./s.
pam.governor_rate_limit = 2;
