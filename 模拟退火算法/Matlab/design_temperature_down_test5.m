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
max_axial_height_R1 = 0.05; % 第一个线圈的最大轴向高度 (m)
max_axial_height_R2 = 0.001; % 第二个线圈的最大轴向高度 (m)
min_r1 = 0.0005; % 第一个线圈的最小导线直径 (m)
max_r1 = 0.002; % 第一个线圈的最大导线直径 (m)
min_r2 = 0.0002; % 第二个线圈的最小导线直径 (m)
max_r2 = 0.001; % 第二个线圈的最大导线直径 (m)
min_n1 = 50; % 第一个线圈的最小环路数量
max_n1 = 2000; % 第一个线圈的最大环路数量
min_n2 = 30; % 第二个线圈的最小环路数量
max_n2 = 1000; % 第二个线圈的最大环路数量


% 模拟退火算法参数
max_iterations = 10000;
initial_temperature = 100000;
final_temperature = 1e-5;
cooling_rate = 0.90;


% 初始化最优参数和最优耦合系数
best_k = 0;
best_R1_distribution = [];
best_R2_distribution = [];
best_d12 = 0;
best_r1 = r1;
best_r2 = r2;
best_n1 = 0;
best_n2 = 0;

% 初始化设计参数
%R1 = max_radius_R1;
%R2 = max_radius_R2;
d12 = 0.12;

% 模拟退火优化过程
temperature = initial_temperature;
iteration = 0;
while temperature > final_temperature && iteration < max_iterations
    
    % 随机生成线圈匝数和导线直径
    n1 = randi([min_n1, max_n1]);
    n2 = randi([min_n2, max_n2]);
    r1 = min_r1 + (max_r1 - min_r1) * rand;
    r2 = min_r2 + (max_r2 - min_r2) * rand;
    
    % 计算环路电感矩阵 L1 和 L2
    L1_matrix = zeros(n1, n1);
    L2_matrix = zeros(n2, n2);
    M_matrix = zeros(n1, n2);
    
    R1_distribution = zeros(n1, 2); % 每匝的半径和轴向高度
    R2_distribution = zeros(n2, 2); % 每匝的半径和轴向高度
    
    current_layer_R1 = min_radius_R1; % 初始化第一层半径（第一个线圈）
    current_height_R1 = 0; % 当前高度（第一个线圈）
    current_layer_R2 = min_radius_R2; % 初始化第一层半径（第二个线圈）
    current_height_R2 = 0; % 当前高度（第二个线圈）
    
    
     for i = 1:n1
        % 确保下一层的每个环路下方有上一层的支持，并且总高度不超过最大轴向高度,且环路不重叠
        if current_height_R1 >= max_axial_height_R1
            break; % 如果高度超过限制，停止增加匝数
        end
        R1_i = current_layer_R1;
        H1_i = current_height_R1; % 当前高度
        R1_distribution(i, :) = [R1_i, H1_i];
        current_height_R1 = current_height_R1 + min_spacing_R1 + r1; % 增加高度，确保环路不重叠
        if current_height_R1 >= max_axial_height_R1 % 如果达到最大高度，增加半径并重置高度
            current_layer_R1 = current_layer_R1 + min_spacing_R1 + r1;
            current_height_R1 = 0;
        end
        for j = 1:n1
            R1_j = R1_distribution(j, 1);
            d_ij = abs(H1_i - R1_distribution(j, 2)); % 环路之间的距离（轴向高度差）
           if R1_i >= min_radius_R1 && R1_i <= max_radius_R1 && d_ij <= max_axial_height_R1
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
        % 确保下一层的每个环路下方有上一层的支持，并且总高度不超过最大轴向高度,且环路不重叠
        if current_height_R2 >= max_axial_height_R2
            break; % 如果高度超过限制，停止增加匝数
        end
        R2_i = current_layer_R2;
        H2_i = current_height_R2; % 当前高度
        R2_distribution(i, :) = [R2_i, H2_i];
        current_height_R2 = current_height_R2 + min_spacing_R2 + r2; % 增加高度，确保环路不重叠
        if current_height_R2 >= max_axial_height_R2 % 如果达到最大高度，增加半径并重置高度
            current_layer_R2 = current_layer_R2 + min_spacing_R2 + r2;
            current_height_R2 = 0;
        end
        for j = 1:n2
            R2_j = R2_distribution(j, 1);
            d_ij = abs(H2_i - R2_distribution(j, 2)); % 环路之间的距离（轴向高度差）
            if R2_i >= min_radius_R2 && R2_i <= max_radius_R2 && d_ij <= max_axial_height_R2
                if i == j
                    L2_matrix(i, j) = mu_0 * R2_i * (log(8 * R2_i / min_spacing_R2) - 7 / 4); % 自感
                else
                    kappa = sqrt((4 * R2_i * R2_j) / ((R2_i + R2_j)^2 + d_ij^2));
                    kappa = min(max(kappa, eps), 1); % 确保 kappa 在 [eps, 1] 范围内，避免 kappa^2 = 0 导致错误
                    [K_kappa, E_kappa] = ellipke(kappa^2); % 完全椭圆积分 K(κ) 和 E(κ)
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
    
     % 计算耦合系数 k，增加惩罚机制，防止匝数过小
    if L1_total > 0 && L2_total > 0
        %penalty = exp(-abs(n1 - 50) / 10) * exp(-abs(n2 - 50) / 10); % 惩罚项，防止匝数过小
        %k = (M_total / sqrt(L1_total * L2_total)) * penalty;
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
    d12_new = d12 + (rand - 0.5) * 0.01 * temperature / initial_temperature;
    d12 = max(0.01, min(0.2, d12_new)); % 确保线圈间距在合理范围内
    
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
fprintf('Best k: %.6f\n', best_k);


% 计算两个线圈的高度、线圈圆环的厚度和最后的线圈半径
height_R1 = max(best_R1_distribution(:, 2)) - min(best_R1_distribution(:, 2)) + min_spacing_R1;
height_R2 = max(best_R2_distribution(:, 2)) - min(best_R2_distribution(:, 2)) + min_spacing_R2;
thickness_R1 = max(best_R1_distribution(:, 1)) - min(best_R1_distribution(:, 1));
thickness_R2 = max(best_R2_distribution(:, 1)) - min(best_R2_distribution(:, 1));
final_radius_R1 = max(best_R1_distribution(:, 1));
final_radius_R2 = max(best_R2_distribution(:, 1));

fprintf('第一个线圈的高度: %.6f m\n', height_R1);
fprintf('第二个线圈的高度: %.6f m\n', height_R2);
fprintf('第一个线圈的厚度: %.6f m\n', thickness_R1);
fprintf('第二个线圈的厚度: %.6f m\n', thickness_R2);
fprintf('第一个线圈的最终半径: %.6f m\n', final_radius_R1);
fprintf('第二个线圈的最终半径: %.6f m\n', final_radius_R2);

% 绘制线圈设计图
figure;
hold on;
axis equal; % 保持图形的纵横比例一致，确保线圈看起来是圆形的

% 绘制第一个线圈 (红色，顶部位置)
for i = 1:best_n1
    theta = linspace(0, 2 * pi, 100);
    x_center = best_R1_distribution(i, 1) * 1000; % 中心位置
    y_center = best_R1_distribution(i, 2) * 1000 + 150; % 高度位置转换为毫米，并调整到图中合适位置
    x = x_center + (r1 / 2) * 1000 * cos(theta); % 半径乘以1000转换为毫米
    y = y_center + (r1 / 2) * 1000 * sin(theta);
    fill(x, y, 'r', 'EdgeColor', 'k');
end

% 绘制第二个线圈 (蓝色，底部位置，对称绘制)
for i = 1:best_n2
    theta = linspace(0, 2 * pi, 100);
    x_center = best_R2_distribution(i, 1) * 1000; % 中心位置
    y_center = best_R2_distribution(i, 2) * 1000 - 150; % 高度位置转换为毫米，并调整到图中合适位置
    x = x_center + (r2 / 2) * 1000 * cos(theta); % 半径乘以1000转换为毫米
    y = y_center + (r2 / 2) * 1000 * sin(theta);
    fill(x, y, 'b', 'EdgeColor', 'k');
end

% 添加线圈之间的连接线
line([0, 0], [-150, 150], 'Color', 'k', 'LineStyle', '-.', 'LineWidth', 1.5);

% 设置图形
xlabel('X (mm)');
ylabel('Y (mm)');
title('线圈设计图');
grid on;
hold off;

