function L = calculateInductance(N, r_m, l_m, d_m)
    % 将输入参数从米转换为英寸
    r = r_m / 0.0254; % 1 inch = 0.0254 meters
    l = l_m / 0.0254;
    d = d_m / 0.0254;

    % 判断是否为单层或空芯线圈
    if d == 0 % 假设当d为0时，认为是单层线圈
        L = (r^2 * N^2) / (8*r + 11*d); % 单层螺旋线圈
    elseif l == 0 % 当l为0时，认为是单层空芯线圈
        L = (r^2 * N^2) / (9*r + 10*d); % 单层空芯线圈
    else % 其他情况认为是多层空芯线圈
        L = (0.8 * r^2 * N^2) / (6*r + 9*l + 10*d); % 多层空芯线圈
    end
end