function [d_guess, eps_guess] = estimate_init_from_phase(fHz, T_meas, d_bounds)
% 基于相位斜率估计 (d*sqrt(eps_r))，给厚度一个合理初值
% d_bounds = [dmin, dmax] (m)
    fHz = fHz(:); T_meas = T_meas(:);
    ph  = unwrap(angle(T_meas));
    % 线性拟合相位 ~ a*f + b
    f = fHz; f = f - mean(f);
    a = (f'*(ph-mean(ph))) / (f'*f + eps);  % rad/Hz  (负号常见)
    % q ≈ (2π f / c0) d sqrt(eps_r)  → d*sqrt(eps_r) ≈ -a*c0
    c0 = 299792458;
    d_sqrt_eps = -a*c0;                   % 单位: m
    % 先给 eps 一个中性值（玻璃/塑料不同，你也可根据文件名自定义）
    eps_guess  = 3.5;                     % 可改：若玻璃可设 5~7
    d_guess    = min(max(d_sqrt_eps / sqrt(eps_guess), d_bounds(1)), d_bounds(2));
    if ~isfinite(d_guess) || d_guess<=0
        d_guess = mean(d_bounds);
    end
end
