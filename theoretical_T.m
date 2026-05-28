function T = theoretical_T(fHz, eps_real, c_param, d_exp, thick_m, theta_deg, pol)
    % ---- 统一为列向量 ----
    fHz = fHz(:);

    eps0 = 8.854187817e-12; mu0 = 4*pi*1e-7; c0 = 1/sqrt(eps0*mu0);
    w = 2*pi*fHz;
    fGHz = fHz/1e9;

    sigma = max(0, c_param) .* max(fGHz,1e-6).^max(d_exp,0);
    eps_c = eps_real - 1j * (sigma ./ (eps0 .* w));

    th   = theta_deg*pi/180;
    sin2 = (sin(th)).^2;

    rt = sqrt(eps_c - sin2);
    eps_c = eps_c + 1e-12*1j;

    q  = (2*pi.*fHz./c0) .* thick_m .* rt;

    switch upper(pol)
        case 'TE'
            denom = (cos(th) + rt);
            Rp = (cos(th) - rt) ./ max(denom, 1e-12);
        case 'TM'
            denom = (eps_c.*cos(th) + rt);
            Rp = (eps_c.*cos(th) - rt) ./ max(denom, 1e-12);
        otherwise
            error('pol must be TE or TM');
    end

    num  = (1 - Rp.^2) .* exp(-1j*q);
    den  = 1 - (Rp.^2) .* exp(-1j*2*q);
    den  = den + 1e-18;

    T = num ./ den;
end
