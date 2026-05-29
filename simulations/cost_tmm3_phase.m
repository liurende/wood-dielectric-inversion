function J = cost_tmm3_phase(x, f, mag_meas_dB, phase_meas, ...
        delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, c0, Z0)
    % Phase-enriched cost function for DE optimization.
    % x = [eps_inf_E, eps_inf_L, f_relax]
    % J = 0.6 * Huber(|S21| error) + 0.4 * wrapped phase MSE

    S21 = tmm_debye_S21(x(1), x(2), x(3), ...
        delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, f, c0, Z0);

    % Magnitude error (Huber loss)
    mag_err = 20 * log10(max(abs(S21), 1e-12)) - mag_meas_dB;
    L_mag = mag_err.^2;
    mask = abs(mag_err) > 1.0;
    L_mag(mask) = 2.0 * abs(mag_err(mask)) - 1.0;

    % Phase error (wrapped to [-pi, pi])
    phase_model = angle(S21);
    phase_err = atan2(sin(phase_model - phase_meas), ...
                      cos(phase_model - phase_meas));
    L_phase = phase_err.^2;

    J = 0.6 * mean(L_mag) + 0.4 * mean(L_phase);
end
