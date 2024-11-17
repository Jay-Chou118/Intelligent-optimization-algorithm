%%通过验证，认为自感矩阵的计算存在问题

%%
N = 11;
r = 150;
l = 30;

% 转换为英寸（1英寸 = 25.4毫米）
mm_to_inch = 25.4;

l_inch = l / mm_to_inch;  % 线圈长度 (inches)
r_inch = r / mm_to_inch;  % 半径 (inches)



% 代入公式计算电感 L (单位：微亨利 μH)
L = ( r_inch^2 * N^2) / (9 * r_inch + 10 * l_inch );
%L = 64.9606


%%
R1 = 150;
R2 = 150;
p = 0;
d = 3;

L1 = totalInductance(11,150,150,3,0)
L1 
%使用遗传算法进行线圈的设计



%%
function M12 = mutualInductance(R1, R2, d12, p)
    % Function to compute the mutual inductance between two coaxial coils
    % using the given formula involving a Bessel function integral.
    %
    % Parameters:
    % R1  - Radius of the first coil (in meters)
    % R2  - Radius of the second coil (in meters)
    % d12 - Distance between the two coils (in meters)
    % p   - A parameter related to the geometry (in meters)
    %
    % Returns:
    % M12 - The mutual inductance (in Henry)

    % Physical constant
    mu0 = 4 * pi * 1e-7; % Permeability of free space (H/m)
    
    % Function to integrate
    integrand = @(x) besselj(1, x * sqrt(R1 / R2)) .* besselj(1, x * sqrt(R2 / R1)) .* besselj(0, x * p / sqrt(R1 * R2)) ...
                    .* exp(-x * d12 / sqrt(R1 * R2));
    
    % Numerical integration from 0 to infinity
    integral_result = integral(integrand, 0, Inf);
    
    % Compute mutual inductance
    M12 = mu0 * pi * sqrt(R1 * R2) * integral_result;
end

function L = selfInductance(R, r)
    % Function to compute the self-inductance of a circular loop
    %
    % Parameters:
    % R - Radius of the loop (in meters)
    % r - Radius of the wire (in meters)
    %
    % Returns:
    % L - The self-inductance (in Henry)

    % Physical constant
    mu0 = 4 * pi * 1e-7; % Permeability of free space (H/m)
    
    % Compute self-inductance
    L = mu0 * R * (log(8 * R / r) - 7 / 4);
end

function L1 = totalInductance(n, R1, R2, d, p)
    % Function to compute the total inductance for n coils
    %
    % Parameters:
    % n  - Number of coils
    % R1 - Radius of the first coil (in meters)
    % R2 - Radius of the second coil (in meters)
    % d  - Distance between the coils (in meters)
    % p  - A parameter related to the geometry (in meters)
    %
    % Returns:
    % L1 - The total inductance (in Henry)

    % Initialize the inductance matrix
    L_matrix = zeros(n, n);
    
    % Fill the inductance matrix
    for i = 1:n
        for j = 1:n
            if i == j
                L_matrix(i, j) = selfInductance(R1, d);
            else
                L_matrix(i, j) = mutualInductance(R1, R2, d, p);
            end
        end
    end
    
    % Calculate the total inductance L1
    L1 = 0;
    for i = 1:n
        L1 = L1 + L_matrix(i, i) + sum(L_matrix(i, [1:i-1, i+1:end]));
    end
end
