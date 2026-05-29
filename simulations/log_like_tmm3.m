function ll = log_like_tmm3(x, f, mag_meas_dB, phase_meas, ...
        delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, c0, Z0)
    % Log-likelihood for MCMC.
    % Proportional to negative phase-enriched cost.
    % For cost-landscape exploration (not calibrated Bayesian inference).
    % Posterior widths reflect relative sensitivity, not physical measurement noise.

    J = cost_tmm3_phase(x, f, mag_meas_dB, phase_meas, ...
        delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, c0, Z0);
    ll = -J;
end
