% 定义常量
mu_0 = 4 * pi * 1e-7; % 真空磁导率 (H/m)

% 定义线圈的参数
n1 = 50; % 第一个线圈的环路数量
n2 = 50; % 第二个线圈的环路数量

% 定义限制条件
min_spacing = 0.005; % 最小环间距 (m)
min_radius_R1 = 0.005; % 第一个线圈的最小半径 (m)
max_radius_R1 = 0.06; % 第一个线圈的最大半径 (m)
min_radius_R2 = 0.002; % 第二个线圈的最小半径 (m)
max_radius_R2 = 0.005; % 第二个线圈的最大半径 (m)
max_axial_height = 0.05; % 最大轴向高度 (m)

% 模拟退火算法参数
max_iterations = 10000;
initial_temperature = 100;
final_temperature = 1e-3;
cooling_rate = 0.95;

% 初始化最优参数和最优耦合系数
best_k = 0;
best_R1 = 0;
best_R2 = 0;
best_d12 = 0;

% 随机初始化设计参数
R1 = min_radius_R1 + (max_radius_R1 - min_radius_R1) * rand;
R2 = min_radius_R2 + (max_radius_R2 - min_radius_R2) * rand;
d12 = 0.12;

% 模拟退火优化过程
temperature = initial_temperature;
iteration = 0;
while temperature > final_temperature && iteration < max_iterations
    % 计算环路电感矩阵 L1 和 L2
    L1_matrix = zeros(n1, n1);
    L2_matrix = zeros(n2, n2);
    M_matrix = zeros(n1, n2);
    
    for i = 1:n1
        R1_i = max_radius_R1 - (i - 1) * min_spacing; % 每个环路的半径根据环间距增加
        for j = 1:n1
            R1_j = R1 + (j - 1) * min_spacing;
            d_ij = abs(i - j) * min_spacing; % 环路之间的距离
            if R1_i >= min_radius_R1 && R1_i <= max_radius_R1 && d_ij <= max_axial_height
                if i == j
                    L1_matrix(i, j) = mu_0 * R1_i * (log(8 * R1_i / min_spacing) - 7 / 4); % 自感
                else
                    kappa = sqrt((4 * R1_i * R1_j) / ((R1_i + R1_j)^2 + d_ij^2));
                    if kappa > 1
                        kappa = 1; % 确保 kappa 在有效范围内
                    end
                    [K_kappa, E_kappa] = ellipke(kappa^2); % 完全椭圆积分 K(κ) 和 E(κ)
                    L1_matrix(i, j) = mu_0 * sqrt(R1_i * R1_j) * ((2 / kappa - kappa) * K_kappa - 2 / kappa * E_kappa); % 互感
                end
            end
        end
    end

    for i = 1:n2
        R2_i = max_radius_R2 - (i - 1) * min_spacing; % 每个环路的半径根据环间距增加
        for j = 1:n2
            R2_j = R2 + (j - 1) * min_spacing;
            d_ij = abs(i - j) * min_spacing; % 环路之间的距离
            if R2_i >= min_radius_R2 && R2_i <= max_radius_R2 && d_ij <= max_axial_height
                if i == j
                    L2_matrix(i, j) = mu_0 * R2_i * (log(8 * R2_i / min_spacing) - 7 / 4); % 自感
                else
                    kappa = sqrt((4 * R2_i * R2_j) / ((R2_i + R2_j)^2 + d_ij^2));
                    if kappa > 1
                        kappa = 1; % 确保 kappa 在有效范围内
                    end
                    [K_kappa, E_kappa] = ellipke(kappa^2); % 完全椭圆积分 K(κ) 和 E(κ)
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
            if R1_i >= min_radius_R1 && R1_i <= max_radius_R1 && R2_j >= min_radius_R2 && R2_j <= max_radius_R2
                kappa = sqrt((4 * R1_i * R2_j) / ((R1_i + R2_j)^2 + d12^2));
                if kappa > 1
                    kappa = 1; % 确保 kappa 在有效范围内
                end
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
        best_R1 = R1;
        best_R2 = R2;
        best_d12 = 0.12;
    end
    
    % 退火过程更新
    R1_new = min_radius_R1 + (max_radius_R1 - min_radius_R1) * (2 * rand - 1) * temperature / initial_temperature;
    R2_new = min_radius_R2 + (max_radius_R2 - min_radius_R2) * (2 * rand - 1) * temperature / initial_temperature;
    d12_new = d12;

    % 确保新值在有效范围内
    R1_new = max(min_radius_R1, min(max_radius_R1, R1_new));
    R2_new = max(min_radius_R2, min(max_radius_R2, R2_new));
    
    % 接受新解的概率
    delta_k = k - best_k;
    if delta_k > 0 || exp(delta_k / temperature) > rand
        R1 = R1_new;
        R2 = R2_new;
        d12 = d12_new;
    end
    
    % 降温
    temperature = temperature * cooling_rate;
    iteration = iteration + 1;
end

% 显示结果
fprintf('最佳耦合系数 k: %.6f\n', best_k);
fprintf('最佳设计参数: R1 = %.6f m, R2 = %.6f m, d12 = %.6f m\n', best_R1, best_R2, best_d12);
