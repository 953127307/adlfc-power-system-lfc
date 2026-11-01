run("A_S1_data_collection.m");

%% 生成控制器权重

m = zeros(n_1,n_2);
G = cell(n_1,n_2); C = cell(n_1,n_2); Y = cell(n_1,n_2);

for i = 1:n_1
    for j = 1:n_2 
        m(i,j) = size(Psi{i,j}, 1); % 基函数维度
        % 中间矩阵
        G{i,j} = [Z{i,j};Psi{i,j}];
        C{i,j} = [1;zeros(m(i,j)-1,1)];
        Y{i,j} = [zeros(2,m(i,j)-1);eye(m(i,j)-1)];
    end
end

k = size(U{1}, 2); % 采样个数

% 检查数据是否满足条件
if ~fun_checkRowFullRank(G)
    error("命题3条件1不满足"); 
end

Psi_ker = cell(n_1,n_2);
for i = 1:n_1
    for j = 1:n_2 
        Psi_ker{i,j} = eye(k)-lsqminnorm(Psi{i,j},Psi{i,j});
        if norm(Z{i,j}*Psi_ker{i,j}) <= sqrt(k)*bar_d(1)
            fprintf("%f, 第%d个区域的第%d个子系统不满足命题3条件2\n", (norm(Z{i,j}*Psi_ker{i,j}))/(sqrt(k)*bar_d(1)),i,j);
            error("0"); 
       end
    end
end

%% 生成权重
p_particular = cell(n_1,n_2);
p = cell(n_1,n_2); Q = cell(n_1,n_2); w = cell(n_1,n_2); 
hat_rho = cell(n_1,n_2); lambda = zeros(n_1,n_2);
for i = 1:n_1
    for j = 1:n_2 
        p_particular{i,j} = lsqminnorm(Psi{i,j},C{i,j});

        karpi_(i,j) = (Z{i,j}*p_particular{i,j}+sqrt(k)*bar_d(1)*norm(p_particular{i,j}))/(norm(Z{i,j}*Psi_ker{i,j})-sqrt(k)*bar_d(1));

        if isfield(pam,'feedback_gain')
            % Algorithm 1 null-space solution, selecting the original gain
            % w(i,j,1)=-1 using recorded next-state/input data only.
            if j < n_2, next_data = X{i,j+1}; else, next_data = U{i}; end
            direction = (Z{i,j}*Psi_ker{i,j})'/norm(Z{i,j}*Psi_ker{i,j});
            karpi(i,j) = (next_data*p_particular{i,j}+pam.feedback_gain(i,j))/(next_data*direction);
            assert(karpi(i,j)>karpi_(i,j),'Selected gain fails robust feasibility');
        end
        p{i,j} = p_particular{i,j} - karpi(i,j)*(Z{i,j}*Psi_ker{i,j})'/norm(Z{i,j}*Psi_ker{i,j});
        row_scale = vecnorm(G{i,j},2,2);
        Q{i,j} = lsqminnorm(G{i,j}./row_scale,Y{i,j}./row_scale);
        if j < n_2
            w{i,j} = (X{i,j+1} * [p{i,j}, Q{i,j}])';
        else
            w{i,j} = (U{i} * [p{i,j}, Q{i,j}])';
        end
        hat_rho{i,j} = -eta(i,j)*w{i,j}(2:end);
        lambda(i,j) = (Z{i,j}-D{i,j})*p{i,j};
        M{i,j} = (Z{i,j}-D{i,j})*Q{i,j};
    end
end
pam.w_11 = w{1,1};
pam.w_12 = w{1,2};
pam.w_13 = w{1,3};
pam.w_14 = w{1,4};

pam.w_21 = w{2,1};
pam.w_22 = w{2,2};
pam.w_23 = w{2,3};
pam.w_24 = w{2,4};

pam.w_31 = w{3,1};
pam.w_32 = w{3,2};
pam.w_33 = w{3,3};
pam.w_34 = w{3,4};

pam.w_41 = w{4,1};
pam.w_42 = w{4,2};
pam.w_43 = w{4,3};
pam.w_44 = w{4,4};

pam.hat_rho_12 = hat_rho{1,2};
pam.hat_rho_13 = hat_rho{1,3};
pam.hat_rho_14 = hat_rho{1,4};

pam.hat_rho_22 = hat_rho{2,2};
pam.hat_rho_23 = hat_rho{2,3};
pam.hat_rho_24 = hat_rho{2,4};

pam.hat_rho_32 = hat_rho{3,2};
pam.hat_rho_33 = hat_rho{3,3};
pam.hat_rho_34 = hat_rho{3,4};

pam.hat_rho_42 = hat_rho{4,2};
pam.hat_rho_43 = hat_rho{4,3};
pam.hat_rho_44 = hat_rho{4,4};

pam.M_12 = M{1,2};
pam.M_13 = M{1,3};
pam.M_14 = M{1,4};

pam.eta = eta;
pam.lambda = lambda;

hat_M_0 = cell(n_1,n_2);
for i = 1:n_1
    for j = 2:n_2
        hat_M_0{i,j} = zeros(m(i,j)-1,1);
    end
end

hat_alpha_0(1,1) = w{1,1}'*[x_0(1,1)-x_d(1)];
hat_alpha_0(1,2) = w{1,2}'*[x_0(1,2)-hat_alpha_0(1,1);sin(x_0(1,1)-x_0(2,1));x_0(1,2);0;1];
hat_alpha_0(1,3) = w{1,3}'*[x_0(1,3)-hat_alpha_0(1,2);x_0(1,3);0];

hat_alpha_0(2,1) = w{2,1}'*[x_0(2,1)-x_d(2)];
hat_alpha_0(2,2) = w{2,2}'*[x_0(2,2)-hat_alpha_0(2,1);sin(x_0(1,1)-x_0(2,1));sin(x_0(2,1)-x_0(3,1));x_0(2,2);0;1];
hat_alpha_0(2,3) = w{2,3}'*[x_0(2,3)-hat_alpha_0(2,2);x_0(2,3);0];

hat_alpha_0(3,1) = w{3,1}'*[x_0(3,1)-x_d(3)];
hat_alpha_0(3,2) = w{3,2}'*[x_0(3,2)-hat_alpha_0(3,1);sin(x_0(2,1)-x_0(3,1));sin(x_0(3,1)-x_0(4,1));x_0(3,2);0;1];
hat_alpha_0(3,3) = w{3,3}'*[x_0(3,3)-hat_alpha_0(3,2);x_0(3,3);0];

hat_alpha_0(4,1) = w{4,1}'*[x_0(4,1)-x_d(4)];
hat_alpha_0(4,2) = w{4,2}'*[x_0(4,2)-hat_alpha_0(4,1);sin(x_0(3,1)-x_0(4,1));x_0(4,2);0;1];
hat_alpha_0(4,3) = w{4,3}'*[x_0(4,3)-hat_alpha_0(4,2);x_0(4,3);0];

pam.hat_alpha_0 = hat_alpha_0;

karpi_(1:4,:)
lambda(1:4,:)

