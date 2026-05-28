function T = theoretical_T_sigma(fHz, eps_real, sigma_f, thick_m, theta_deg, pol)
% 与 theoretical_T 相同公式，但损耗项直接传入 sigma(f) 向量
    fHz = fHz(:);
    sigma_f = sigma_f(:);
    eps0 = 8.854187817e-12; mu0 = 4*pi*1e-7; c0 = 1/sqrt(eps0*mu0);
    w = 2*pi*fHz;

    eps_c = eps_real - 1j*(sigma_f ./ (eps0 .* w + 1e-30));  % 防0
    th   = theta_deg*pi/180;  sin2 = (sin(th)).^2;
    rt   = sqrt(eps_c - sin2);
    q    = (2*pi.*fHz./c0) .* thick_m .* rt;

    switch upper(pol)
        case 'TE'
            Rp = (cos(th) - rt) ./ (cos(th) + rt + 1e-18);
        case 'TM'
            Rp = (eps_c.*cos(th) - rt) ./ (eps_c.*cos(th) + rt + 1e-18);
        otherwise
            error('pol must be TE or TM');
    end

    num = (1 - Rp.^2) .* exp(-1j*q);
    den = 1 - (Rp.^2) .* exp(-1j*2*q) + 1e-18;
    T = num ./ den;
end
