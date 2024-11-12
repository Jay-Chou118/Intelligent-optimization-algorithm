clear; clc;

% 定义常量
mu_0 = 4 * pi * 1e-7; % 真空磁导率 (H/m)
rho = 1.68e-8; % 铜的电阻率 (ohm*m)

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
max_radius_R2 = 0.006; % 第二个线圈的最大半径 (m)
max_axial_height_R1 = 0.05; % 第一个线圈的最大轴向高度 (m)
max_axial_height_R2 = 0.01; % 第二个线圈的最大轴向高度 (m)
min_r1 = 0.0005; % 第一个线圈的最小导线直径 (m)
max_r1 = 0.002; % 第一个线圈的最大导线直径 (m)
min_r2 = 0.0002; % 第二个线圈的最小导线直径 (m)
max_r2 = 0.001; % 第二个线圈的最大导线直径 (m)
min_n1 = 50; % 第一个线圈的最小环路数量
max_n1 = 2000; % 第一个线圈的最大环路数量
min_n2 = 30; % 第二个线圈的最小环路数量
max_n2 = 1000; % 第二个线圈的最大环路数量

% 设置多目标遗传算法 (MOGA) 参数
population_size = 100; % 种群大小
max_generations = 200; % 最大迭代代数

% 定义优化问题
objective_function = @(x) coil_objectives(x, mu_0, rho, min_spacing_R1, min_spacing_R2);

nvars = 4; % 决策变量的数量 (n1, n2, r1, r2)
lb = [min_n1, min_n2, min_r1, min_r2]; % 决策变量的下界
ub = [max_n1, max_n2, max_r1, max_r2]; % 决策变量的上界

% 运行多目标遗传算法 (MOGA)
options = optimoptions('ga', 'PopulationSize', population_size, 'MaxGenerations', max_generations, 'Display', 'iter', 'PlotFcn', {@gaplotpareto}, 'UseParallel', true);
[x_opt, fval] = ga(objective_function, nvars, [], [], [], [], lb, ub, [], options);

% 显示结果
fprintf('优化后的设计参数:\n');
fprintf('最佳第一个线圈环路数量 n1 = %d\n', round(x_opt(1)));
fprintf('最佳第二个线圈环路数量 n2 = %d\n', round(x_opt(2)));
fprintf('最佳第一个线圈导线直径 r1 = %.6f m\n', x_opt(3));
fprintf('最佳第二个线圈导线直径 r2 = %.6f m\n', x_opt(4));

% 定义目标函数
function f = coil_objectives(x, mu_0, rho, min_spacing_R1, min_spacing_R2)
    % 提取设计变量
    n1 = round(x(1));
    n2 = round(x(2));
    r1 = x(3);
    r2 = x(4);

    % 计算线圈匝的分布和几何参数
    min_radius_R1 = 0.005;
    max_radius_R1 = 0.06;
    min_radius_R2 = 0.002;
    max_radius_R2 = 0.006;
    max_axial_height_R1 = 0.05;
    max_axial_height_R2 = 0.01;
    
    % 第一个线圈的几何参数计算
    R1_distribution = calculate_coil_distribution(n1, r1, min_radius_R1, max_axial_height_R1, min_spacing_R1);
    R2_distribution = calculate_coil_distribution(n2, r2, min_radius_R2, max_axial_height_R2, min_spacing_R2);

    % 计算自感和互感
    L1_total = 0;
    L2_total = 0;
    M_total = 0;
    d12 = 0.12;
    for i = 1:n1
        R1_i = R1_distribution(i, 1);
        L1_total = L1_total + (mu_0 * R1_i * n1^2) / (8 * R1_i + 11 * min_spacing_R1);
        for j = 1:n2
            R2_j = R2_distribution(j, 1);
            kappa = sqrt((4 * R1_i * R2_j) / ((R1_i + R2_j)^2 + d12^2));
            [K_kappa, E_kappa] = ellipke(kappa^2);
            M_total = M_total + mu_0 * sqrt(R1_i * R2_j) * ((2 / kappa - kappa) * K_kappa - 2 / kappa * E_kappa);
        end
    end

    for i = 1:n2
        R2_i = R2_distribution(i, 1);
        L2_total = L2_total + (mu_0 * R2_i * n2^2) / (8 * R2_i + 11 * min_spacing_R2);
    end

    % 计算电阻 R1 和 R2
    length_R1 = 2 * pi * mean(R1_distribution(:, 1)) * n1;
    length_R2 = 2 * pi * mean(R2_distribution(:, 1)) * n2;
    A1 = pi * (r1 / 2)^2;
    A2 = pi * (r2 / 2)^2;
    R1 = rho * length_R1 / A1;
    R2 = rho * length_R2 / A2;

    % 计算耦合系数 k
    if L1_total > 0 && L2_total > 0
        k = M_total / sqrt(L1_total * L2_total);
    else
        k = 0;
    end

    % 定义多目标优化函数：最大化 -k，最小化 R1 和 R2
    f = [-k, R1, R2];
end

% 计算线圈分布的辅助函数
function distribution = calculate_coil_distribution(n, r, min_radius, max_height, min_spacing)
    distribution = zeros(n, 2);
    current_layer = min_radius;
    current_height = 0;
    for i = 1:n
        distribution(i, :) = [current_layer, current_height];
        current_height = current_height + min_spacing + r;
        if current_height >= max_height
            current_layer = current_layer + min_spacing + r;
            current_height = 0;
        end
    end
end
