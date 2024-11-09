clear;clc;
% 定义常量
mu_0 = 4 * pi * 1e-7; % 真空磁导率 (H/m)

% 定义线圈的参数
R1 = 0.1; % 第一个线圈的初始半径 (m)
R2 = 0.15; % 第二个线圈的初始半径 (m)
r1 = 0.01; % 第一个线圈的导线半径 (m)
r2 = 0.015; % 第二个线圈的导线半径 (m)
d12 = 0.12; % 两个线圈之间的距离 (m)
n1 = 10; % 第一个线圈的环路数量
n2 = 10; % 第二个线圈的环路数量


% 定义限制条件
min_spacing = 0.005; % 最小环间距 (m)
min_radius = 0.05; % 最小线圈半径 (m)
max_radius = 0.2; % 最大线圈半径 (m)
max_axial_height = 0.3; % 最大轴向高度 (m)


% 初始化环路电感矩阵 L1 和 L2
L1_matrix = zeros(n1, n1);
L2_matrix = zeros(n2, n2);
M_matrix = zeros(n1, n2);

% 计算第一个线圈的环路电感矩阵 L1
for i = 1:n1
    R1_i = R1 + (i - 1) * min_spacing; % 每个环路的半径根据环间距增加
    for j = 1:n1
        R1_j = R1 + (j - 1) * min_spacing;
        d_ij = abs(i - j) * min_spacing; % 环路之间的距离
        if R1_i >= min_radius && R1_i <= max_radius && d_ij <= max_axial_height
            if i == j
                L1_matrix(i, j) = mu_0 * R1_i * (log(8 * R1_i / r1) - 7 / 4); % 自感
            else
                kappa = sqrt((4 * R1_i * R1_j) / ((R1_i + R1_j)^2 + d_ij^2));
                K_kappa = ellipke(kappa^2); % 完全椭圆积分 K(κ)
                E_kappa = ellipke(1 - kappa^2); % 完全椭圆积分 E(κ)
                L1_matrix(i, j) = mu_0 * sqrt(R1_i * R1_j) * ((2 / kappa - kappa) * K_kappa - 2 / kappa * E_kappa); % 互感
            end
        end
    end
end

% 计算第二个线圈的环路电感矩阵 L2
for i = 1:n2
    R2_i = R2 + (i - 1) * min_spacing; % 每个环路的半径根据环间距增加
    for j = 1:n2
        R2_j = R2 + (j - 1) * min_spacing;
        d_ij = abs(i - j) * min_spacing; % 环路之间的距离
        if R2_i >= min_radius && R2_i <= max_radius && d_ij <= max_axial_height
            if i == j
                L2_matrix(i, j) = mu_0 * R2_i * (log(8 * R2_i / r2) - 7 / 4); % 自感
            else
                kappa = sqrt((4 * R2_i * R2_j) / ((R2_i + R2_j)^2 + d_ij^2));
                K_kappa = ellipke(kappa^2); % 完全椭圆积分 K(κ)
                E_kappa = ellipke(1 - kappa^2); % 完全椭圆积分 E(κ)
                L2_matrix(i, j) = mu_0 * sqrt(R2_i * R2_j) * ((2 / kappa - kappa) * K_kappa - 2 / kappa * E_kappa); % 互感
            end
        end
    end
end

% 计算两个线圈之间的互感矩阵 M
for i = 1:n1
    R1_i = R1 + (i - 1) * min_spacing; % 每个环路的半径根据环间距增加
    for j = 1:n2
        R2_j = R2 + (j - 1) * min_spacing;
        if R1_i >= min_radius && R1_i <= max_radius && R2_j >= min_radius && R2_j <= max_radius && d12 <= max_axial_height
            kappa = sqrt((4 * R1_i * R2_j) / ((R1_i + R2_j)^2 + d12^2));
            K_kappa = ellipke(kappa^2); % 完全椭圆积分 K(κ)
            E_kappa = ellipke(1 - kappa^2); % 完全椭圆积分 E(κ)
            M_matrix(i, j) = mu_0 * sqrt(R1_i * R2_j) * ((2 / kappa - kappa) * K_kappa - 2 / kappa * E_kappa);
        end
    end
end



% 计算第一个线圈的总自感 L1
L1_total = sum(diag(L1_matrix)) + sum(sum(L1_matrix)) - sum(diag(L1_matrix));

% 计算第二个线圈的总自感 L2
L2_total = sum(diag(L2_matrix)) + sum(sum(L2_matrix)) - sum(diag(L2_matrix));

% 计算两个线圈之间的总互感 M
M_total = sum(sum(M_matrix));

% 计算耦合系数 k
k = M_total / sqrt(L1_total * L2_total);

% 显示结果
fprintf('第一个线圈的自感矩阵 L1:\n');
disp(L1_matrix);

fprintf('第二个线圈的自感矩阵 L2:\n');
disp(L2_matrix);

fprintf('两个线圈之间的互感矩阵 M\n');
disp(M_matrix);


fprintf('第一个线圈的总自感 L1_total: %.6f mH\n', L1_total * 1000);
fprintf('第二个线圈的总自感 L2_total: %.6f mH\n', L2_total * 1000);
fprintf('两个线圈之间的总互感 M_total: %.6f mH\n', M_total * 1000);

