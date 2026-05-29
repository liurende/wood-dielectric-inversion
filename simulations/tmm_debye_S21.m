function S21 = tmm_debye_S21(eps_inf_E, eps_inf_L, f_relax, ...
        delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, f, c0, Z0)
    % Compute S21 for alternating E-L layered structure.
    % Each layer uses Debye dispersion.
    % Uses Chebyshev acceleration: O(1) per frequency, independent of N_pairs.

    eps_E = eps_inf_E + delta_eps_E ./ (1 + 1j * f / f_relax);
    eps_L = eps_inf_L + delta_eps_L ./ (1 + 1j * f / f_relax);

    [A, B, C, D] = tmm_chebyshev(eps_E, eps_L, d_E, d_L, N_pairs, f, c0, Z0);

    denom = A + B/Z0 + C*Z0 + D;
    S21 = 2 ./ denom;
end

function [A, B, C, D] = tmm_chebyshev(eps_E, eps_L, d_E, d_L, N_pairs, f, c0, Z0)
    % Compute total ABCD matrix using Chebyshev polynomial identity.
    % All inputs are N_freq x 1 vectors (eps_E, eps_L may be complex).
    % d_E, d_L are scalars.
    % Returns A, B, C, D as N_freq x 1 column vectors.

    n_E = sqrt(eps_E);
    n_L = sqrt(eps_L);

    k0 = 2 * pi * f / c0;

    delta_E = k0 .* n_E * d_E;
    delta_L = k0 .* n_L * d_L;

    Z_E = Z0 ./ n_E;
    Z_L = Z0 ./ n_L;

    cE = cos(delta_E); sE = sin(delta_E);
    cL = cos(delta_L); sL = sin(delta_L);

    alpha = cE.*cL - (Z_E./Z_L).*sE.*sL;
    beta  = 1j * (Z_E.*sE.*cL + Z_L.*cE.*sL);
    gamma = 1j * ((1./Z_E).*sE.*cL + (1./Z_L).*cE.*sL);
    delta_term = cE.*cL - (Z_L./Z_E).*sE.*sL;

    s = (alpha + delta_term) / 2;
    theta = acos(s);

    sin_theta = sin(theta);
    near_degenerate = abs(sin_theta) < 1e-7;

    U_Nm1 = zeros(size(s));
    U_Nm2 = zeros(size(s));

    if any(near_degenerate)
        sgn = sign(real(s(near_degenerate)));
        U_Nm1(near_degenerate) = N_pairs * sgn.^(N_pairs - 1);
        U_Nm2(near_degenerate) = (N_pairs - 1) * sgn.^(N_pairs - 2);
    end

    ok = ~near_degenerate;
    U_Nm1(ok) = sin(N_pairs * theta(ok)) ./ sin_theta(ok);
    U_Nm2(ok) = sin((N_pairs - 1) * theta(ok)) ./ sin_theta(ok);

    A = U_Nm1 .* alpha      - U_Nm2;
    B = U_Nm1 .* beta;
    C = U_Nm1 .* gamma;
    D = U_Nm1 .* delta_term - U_Nm2;
end
