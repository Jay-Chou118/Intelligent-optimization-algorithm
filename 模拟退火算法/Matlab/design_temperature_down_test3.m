clear;clc;
% 定义常量
mu_0 = 4 * pi * 1e-7; % 真空磁导率 (H/m)
% 定义线圈的初始参数
n1 = 50; % 第一个线圈的初始环路数量
n2 = 50; % 第二个线圈的初始环路数量
r1 = 0.001; % 第一个线圈的初始导线直径 (m)
r2 = 0.0008; % 第二个线圈的初始导线直径 (m)

% 定义限制条件
min_spacing_R1 = r1; % 第一个线圈的最小环间距为导线直径 (m)
min_spacing_R2 = r2; % 第二个线圈的最小环间距为导线直径 (m)
min_radius_R1 = 0.005; % 第一个线圈的最小半径 (m)
max_radius_R1 = 0.06; % 第一个线圈的最大半径 (m)
min_radius_R2 = 0.002; % 第二个线圈的最小半径 (m)
max_radius_R2 = 0.005; % 第二个线圈的最大半径 (m)
max_axial_height = 0.05; % 最大轴向高度 (m)
min_r1 = 0.0005; % 第一个线圈的最小导线直径 (m)
max_r1 = 0.002; % 第一个线圈的最大导线直径 (m)
min_r2 = 0.0002; % 第二个线圈的最小导线直径 (m)
max_r2 = 0.001; % 第二个线圈的最大导线直径 (m)
min_n1 = 5; % 第一个线圈的最小环路数量
max_n1 = 100; % 第一个线圈的最大环路数量
min_n2 = 10; % 第二个线圈的最小环路数量
max_n2 = 100; % 第二个线圈的最大环路数量


% 模拟退火算法参数
max_iterations = 10000;
initial_temperature = 100000;
final_temperature = 1e-5;
cooling_rate = 0.999;


% 初始化最优参数和最优耦合系数
best_k = 0;
best_R1_distribution = zeros(n1, 2); % 存储第一个线圈每匝的半径和轴向高度
best_R2_distribution = zeros(n2, 2); % 存储第二个线圈每匝的半径和轴向高度
best_d12 = 0;
best_r1 = r1;
best_r2 = r2;
best_n1 = n1;
best_n2 = n2;

% 初始化设计参数
R1 = max_radius_R1;
R2 = max_radius_R2;
d12 = 0.12;

% 模拟退火优化过程
temperature = initial_temperature;
iteration = 0;
while temperature > final_temperature && iteration < max_iterations
    % 计算环路电感矩阵 L1 和 L2
    L1_matrix = zeros(n1, n1);
    L2_matrix = zeros(n2, n2);
    M_matrix = zeros(n1, n2);
    
    R1_distribution = zeros(n1, 2); % 每匝的半径和轴向高度
    R2_distribution = zeros(n2, 2); % 每匝的半径和轴向高度
    
    for i = 1:n1
        R1_i = max_radius_R1 - (i - 1) * min_spacing_R1; % 每个环路的半径根据环间距递减
        H1_i = (i - 1) * min_spacing_R1; % 轴向高度
        R1_distribution(i, :) = [R1_i, H1_i];
        for j = 1:n1
            R1_j = max_radius_R1 - (j - 1) * min_spacing_R1;
            d_ij = abs(i - j) * min_spacing_R1; % 环路之间的距离
            if R1_i >= min_radius_R1 && R1_i <= max_radius_R1 && d_ij <= max_axial_height
                if i == j
                    L1_matrix(i, j) = mu_0 * R1_i * (log(8 * R1_i / min_spacing_R1) - 7 / 4); % 自感
                else
                    kappa = sqrt((4 * R1_i * R1_j) / ((R1_i + R1_j)^2 + d_ij^2));
                    kappa = min(max(kappa, eps), 1); % 确保 kappa 在 [eps, 1] 范围内，避免 kappa^2 = 0 导致错误
                    [K_kappa, E_kappa] = ellipke(max(0, min(1, kappa^2))); % 确保 kappa^2 在 [0, 1] 范围内，计算完全椭圆积分 K(κ) 和 E(κ)
                    L1_matrix(i, j) = mu_0 * sqrt(R1_i * R1_j) * ((2 / kappa - kappa) * K_kappa - 2 / kappa * E_kappa); % 互感
                end
            end
        end
    end

    for i = 1:n2
        R2_i = max_radius_R2 - (i - 1) * min_spacing_R2; % 每个环路的半径根据环间距递减
        H2_i = (i - 1) * min_spacing_R2; % 轴向高度
        R2_distribution(i, :) = [R2_i, H2_i];
        for j = 1:n2
            R2_j = max_radius_R2 - (j - 1) * min_spacing_R2;
            d_ij = abs(i - j) * min_spacing_R2; % 环路之间的距离
            if R2_i >= min_radius_R2 && R2_i <= max_radius_R2 && d_ij <= max_axial_height
                if i == j
                    L2_matrix(i, j) = mu_0 * R2_i * (log(8 * R2_i / min_spacing_R2) - 7 / 4); % 自感
                else
                    kappa = sqrt((4 * R2_i * R2_j) / ((R2_i + R2_j)^2 + d_ij^2));
                    kappa = min(max(kappa, eps), 1); % 确保 kappa 在 [eps, 1] 范围内，避免 kappa^2 = 0 导致错误
                    [K_kappa, E_kappa] = ellipke(max(0, min(1, kappa^2))); % 确保 kappa^2 在 [0, 1] 范围内，计算完全椭圆积分 K(κ) 和 E(κ)
                    L2_matrix(i, j) = mu_0 * sqrt(R2_i * R2_j) * ((2 / kappa - kappa) * K_kappa - 2 / kappa * E_kappa); % 互感
                end
            end
        end
    end

    % 计算两个线圈之间的互感矩阵 M
    for i = 1:n1
        R1_i = R1_distribution(i, 1); % 使用分布中的半径
        for j = 1:n2
            R2_j = R2_distribution(j, 1); % 使用分布中的半径
            if R1_i >= min_radius_R1 && R1_i <= max_radius_R1 && R2_j >= min_radius_R2 && R2_j <= max_radius_R2
                kappa = sqrt((4 * R1_i * R2_j) / ((R1_i + R2_j)^2 + d12^2));
                kappa = min(max(kappa, eps), 1); % 确保 kappa 在 [eps, 1] 范围内，避免 kappa^2 = 0 导致错误
                [K_kappa, E_kappa] = ellipke(kappa^2); % 完全椭圆积分 K(κ) 和 E(κ)
                M_matrix(i, j) = mu_0 * sqrt(R1_i * R2_j) * ((2 / kappa - kappa) * K_kappa - 2 / kappa * E_kappa);
            end
        end
    end

    % 计算总自感和互感
    L1_total = sum(L1_matrix(:));
    L2_total = sum(L2_matrix(:));
    M_total = sum(M_matrix(:));
    
    % 计算耦合系数 k
    if L1_total > 0 && L2_total > 0
        k = M_total / sqrt(L1_total * L2_total);
    else
        k = 0;
    end
    
    % 更新最优解
    if k > best_k
        best_k = k;
        best_R1_distribution = R1_distribution;
        best_R2_distribution = R2_distribution;
        best_d12 = 0.12;
        best_r1 = r1;
        best_r2 = r2;
        best_n1 = n1;
        best_n2 = n2;
    end
    
    % 退火过程更新
    R1_new = R1 + (rand - 0.5) * (max_radius_R1 - min_radius_R1) * temperature / initial_temperature;
    R2_new = R2 + (rand - 0.5) * (max_radius_R2 - min_radius_R2) * temperature / initial_temperature;
    r1_new = r1 + (rand - 0.5) * (max_r1 - min_r1) * temperature / initial_temperature;
    r2_new = r2 + (rand - 0.5) * (max_r2 - min_r2) * temperature / initial_temperature;
    n1_new = round(n1 + (rand - 0.5) * (max_n1 - min_n1) * temperature / initial_temperature);
    n2_new = round(n2 + (rand - 0.5) * (max_n2 - min_n2) * temperature / initial_temperature);

    % 确保新值在有效范围内
    R1_new = max(min_radius_R1, min(max_radius_R1, R1_new));
    R2_new = max(min_radius_R2, min(max_radius_R2, R2_new));
    r1_new = max(min_r1, min(max_r1, r1_new));
    r2_new = max(min_r2, min(max_r2, r2_new));
    n1_new = max(min_n1, min(max_n1, n1_new));
    n2_new = max(min_n2, min(max_n2, n2_new));
    
    % 接受新解的概率
    delta_k = k - best_k;
    if delta_k > 0 || exp(delta_k / temperature) > rand
        R1 = R1_new;
        R2 = R2_new;
        r1 = r1_new;
        r2 = r2_new;
        n1 = n1_new;
        n2 = n2_new;
    end
    
    % 降温
    temperature = temperature * cooling_rate;
    iteration = iteration + 1;
    fprintf('Iteration %d, Temperature %.4f, Best k: %.6f\n', iteration, temperature, best_k);
end

% 显示结果
fprintf('最佳耦合系数 k: %.6f\n', best_k);
fprintf('最佳设计参数:\n');
fprintf('第一个线圈每匝的半径和轴向高度:\n');
disp(best_R1_distribution);
fprintf('第二个线圈每匝的半径和轴向高度:\n');
disp(best_R2_distribution);
fprintf('最佳线圈之间的距离 d12 = %.6f m\n', best_d12);
fprintf('最佳第一个线圈导线直径 r1 = %.6f m\n', best_r1);
fprintf('最佳第二个线圈导线直径 r2 = %.6f m\n', best_r2);
fprintf('最佳第一个线圈环路数量 n1 = %d\n', best_n1);
fprintf('最佳第二个线圈环路数量 n2 = %d\n', best_n2);
fprintf('第一个线圈的总自感 L1: %.6e mH\n', L1_total*1000);
fprintf('第二个线圈的总自感 L2: %.6e mH\n', L2_total*1000);
fprintf('两个线圈之间的互感 M: %.6e H\n', M_total);
