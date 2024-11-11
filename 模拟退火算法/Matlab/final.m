clear;clc;
% 定义常量
mu_0 = 4 * pi * 1e-7; % 真空磁导率 (H/m)
rho = 1.68e-8; % 铜的电阻率 (ohm*m)
% 定义线圈的初始参数
n1 = 50; % 第一个线圈的初始环路数量
n2 = 50; % 第二个线圈的初始环路数量
r1 = 0.001; % 第一个线圈的初始导线直径 (m)
r2 = 0.0008; % 第二个线圈的初始导线直径 (m)

% 定义限制条件
min_spacing_R1 = r1; % 更新最小匝间距以确保其随导线直径变化 % 第一个线圈的最小环间距为导线直径 (m)
min_spacing_R2 = r2; % 更新最小匝间距以确保其随导线直径变化 % 第二个线圈的最小环间距为导线直径 (m)
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
max_layers_R1 = floor(max_axial_height_R1 / r1); % 第一个线圈的最大层数
max_layers_R2 = floor(max_axial_height_R2 / r2); % 第二个线圈的最大层数

% 模拟退火算法参数
max_iterations = 100000;
initial_temperature = 100000;
final_temperature = 1e-5;
cooling_rate = 0.95;

% 初始化最优参数和最优耦合系数
best_k = 0;
best_R1_distribution = [];
best_R2_distribution = [];
best_d12 = 0;
best_r1 = r1;
best_r2 = r2;
best_n1 = 0;
best_n2 = 0;
best_R = Inf;

% 初始化设计参数
d12 = 0.12;

% 模拟退火优化过程
temperature = initial_temperature;
iteration = 0;
while temperature > final_temperature && iteration < max_iterations
    % 随机生成线圈匝数和导线直径
    step_size_factor = temperature / initial_temperature; % 自适应步长因子，随着温度降低逐渐减小
    n1 = round(min_n1 + (max_n1 - min_n1) * step_size_factor * rand);
    n1 = max(min(n1, max_n1), min_n1);
    n2 = round(min_n2 + (max_n2 - min_n2) * step_size_factor * rand);
    n2 = max(min(n2, max_n2), min_n2);
    r1 = min_r1 + (max_r1 - min_r1) * step_size_factor * rand;
    r1 = max(min(r1, max_r1), min_r1);
    r2 = min_r2 + (max_r2 - min_r2) * step_size_factor * rand;
    r2 = max(min(r2, max_r2), min_r2);
    
    % 计算线圈匝的分布
    R1_distribution = zeros(n1, 2); % 每匝的半径和轴向高度
    R2_distribution = zeros(n2, 2); % 每匝的半径和轴向高度
    current_layer_R1 = min_radius_R1; % 初始化第一层半径（第一个线圈）
    current_height_R1 = 0; % 当前高度（第一个线圈）
    for i = 1:n1
        R1_distribution(i, :) = [current_layer_R1, current_height_R1];
        current_height_R1 = current_height_R1 + min_spacing_R1 + r1;
        if current_height_R1 >= max_axial_height_R1 % 如果达到最大高度，增加半径并重置高度
            current_layer_R1 = current_layer_R1 + min_spacing_R1 + r1;
            current_height_R1 = 0;
        end
    end
    
    current_layer_R2 = min_radius_R2; % 初始化第一层半径（第二个线圈）
    current_height_R2 = 0; % 当前高度（第二个线圈）
    for i = 1:n2
        R2_distribution(i, :) = [current_layer_R2, current_height_R2];
        current_height_R2 = current_height_R2 + min_spacing_R2 + r2;
        if current_height_R2 >= max_axial_height_R2 % 如果达到最大高度，增加半径并重置高度
            current_layer_R2 = current_layer_R2 + min_spacing_R2 + r2;
            current_height_R2 = 0;
        end
    end
    
<<<<<<< HEAD
    % 计算自感和互感
    L1_total = 0;
    L2_total = 0;
    M_total = 0;
=======

    % 使用矢量化和并行计算来优化自感矩阵 L1 的计算
% 使用矢量化和并行计算来优化互感矩阵 M 的计算parfor i = 1:n1
    R1_i = R1_distribution(i, 1);
    R2_distribution_local = R2_distribution; % 创建本地版本，减少对共享变量的依赖
    for j = 1:n2
        R2_j = R2_distribution_local(j, 1);
        if R1_i >= min_radius_R1 && R1_i <= max_radius_R1 && R2_j >= min_radius_R2 && R2_j <= max_radius_R2
            syms phi1 phi2 R1_sym R2_sym d_sym
            R1_sym = R1_i;
            R2_sym = R2_j;
            d_sym = d12;
            M_sym = mu_0 * R1_sym * R2_sym / sqrt((R1_sym * cos(phi1) - R2_sym * cos(phi2))^2 + (R1_sym * sin(phi1) - R2_sym * sin(phi2))^2 + d_sym^2);
            M_matrix(i, j) = double(vpa(int(int(M_sym, phi1, 0, 2 * pi), phi2, 0, 2 * pi)));
        end
    end
end

    % 计算两个线圈之间的互感矩阵 M
>>>>>>> b1203998b8afd0d72bae25d3c01e1f9871170782
    for i = 1:n1
        R1_i = R1_distribution(i, 1);
        % 使用 Grover 方法计算自感
        L1_total = L1_total + (mu_0 * R1_i * n1^2) / (8 * R1_i + 11 * min_spacing_R1);
        for j = 1:n2
            R2_j = R2_distribution(j, 1);
            % 计算互感
            d_ij = d12;
            kappa = sqrt((4 * R1_i * R2_j) / ((R1_i + R2_j)^2 + d_ij^2));
            [K_kappa, E_kappa] = ellipke(kappa^2);
            M_total = M_total + mu_0 * sqrt(R1_i * R2_j) * ((2 / kappa - kappa) * K_kappa - 2 / kappa * E_kappa);
        end
    end
    
    for i = 1:n2
        R2_i = R2_distribution(i, 1);
        % 使用 Grover 方法计算自感
        L2_total = L2_total + (mu_0 * R2_i * n2^2) / (8 * R2_i + 11 * min_spacing_R2);
    end
    
    % 计算电阻 R1 和 R2
    length_R1 = 2 * pi * mean(R1_distribution(:, 1)) * n1; % 第一个线圈导线的总长度
    length_R2 = 2 * pi * mean(R2_distribution(:, 1)) * n2; % 第二个线圈导线的总长度
    A1 = pi * (r1 / 2)^2; % 第一个线圈导线的截面积
    A2 = pi * (r2 / 2)^2; % 第二个线圈导线的截面积
    R1 = rho * length_R1 / A1; % 第一个线圈的电阻
    R2 = rho * length_R2 / A2; % 第二个线圈的电阻
    total_R = R1 + R2; % 总电阻
    
    % 计算耦合系数 k，增加惩罚机制，防止匝数过小
    if L1_total > 0 && L2_total > 0
        % 增加惩罚机制，防止不符合几何约束的解
        penalty_factor = 1;
        if current_height_R1 > max_axial_height_R1 || current_layer_R1 > max_radius_R1 || current_height_R2 > max_axial_height_R2 || current_layer_R2 > max_radius_R2
            penalty_factor = 0.5; % 如果违反几何约束，施加惩罚
        end
        k = (M_total / sqrt(L1_total * L2_total)) * penalty_factor;
    else
        k = 0;
    end
    
    % 多目标优化：考虑最大化耦合系数 k，同时最小化电阻 R
    w_k = 0.5; % 权重，用于平衡耦合系数和电阻的优化
    w_R = 0.5;
    objective_value = w_k * k - w_R * total_R; % 最大化 k，最小化 R
    
    % 更新最优解
    if objective_value > (w_k * best_k - w_R * best_R)
        best_k = k;
        best_R = total_R;
        best_R1_distribution = R1_distribution;
        best_R2_distribution = R2_distribution;
        best_d12 = 0.12;
        best_r1 = r1;
        best_r2 = r2;
        best_n1 = n1;
        best_n2 = n2;
    end
    
    % 降温
    temperature = temperature / (1 + 0.001 * iteration); % 使用自适应冷却策略
    iteration = iteration + 1;
    fprintf('Iteration %d, Temperature %.4f, Best k: %.6f, Best R: %.6f\n', iteration, temperature, best_k, best_R);
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
fprintf('第一个线圈的总自感 L1: %.6e mH\n', L1_total * 1000);
fprintf('第二个线圈的总自感 L2: %.6e mH\n', L2_total * 1000);
fprintf('两个线圈之间的互感 M: %.6e H\n', M_total);
fprintf('总电阻 R: %.6f ohm\n', best_R);
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
% 绘制第一个线圈 (红色，顶部位置)
scatter(-best_R1_distribution(:, 1) * 1000, best_R1_distribution(:, 2) * 1000 + 150, 'r', 'filled');
scatter(best_R1_distribution(:, 1) * 1000, best_R1_distribution(:, 2) * 1000 + 150, 'r', 'filled');
rectangle('Position',[-max(best_R1_distribution(:, 1)) * 1000, 150, 2 * max(best_R1_distribution(:, 1)) * 1000, 50], 'EdgeColor', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
% 绘制第二个线圈 (蓝色，底部位置)
scatter(-best_R2_distribution(:, 1) * 1000, best_R2_distribution(:, 2) * 1000 - 150, 'b', 'filled');
scatter(best_R2_distribution(:, 1) * 1000, best_R2_distribution(:, 2) * 1000 - 150, 'b', 'filled');
rectangle('Position',[-max(best_R2_distribution(:, 1)) * 1000, -200, 2 * max(best_R2_distribution(:, 1)) * 1000, 50], 'EdgeColor', 'b', 'LineStyle', '--', 'LineWidth', 1.5);
% 添加线圈之间的连接线
line([0, 0], [-150, 150], 'Color', 'k', 'LineStyle', '-.', 'LineWidth', 1.5);
% 设置图形属性
xlabel('Lateral direction (mm)');
ylabel('Axial direction (mm)');
title('Detailed Coil Design with Separation');
legend('Coil 1 (Top)', 'Coil 2 (Bottom)', 'Location', 'Best');
grid on;
axis equal;
hold off;

