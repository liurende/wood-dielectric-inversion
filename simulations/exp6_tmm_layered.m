%% Experiment 6: Transfer Matrix Method for Layered Wood Model
% Alternating earlywood/latewood layers to capture annual ring effects.
% Key question: does the layered structure create identifiable spectral
% features that the uniform Fabry-Perot model cannot reproduce?

function exp6_tmm_layered()
    clc; close all;

    outDir = fullfile(pwd, 'simulations', 'results');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    fontName = 'Times New Roman';
    c0 = 299792458;
    Z0 = 377;  % free-space wave impedance

    fGHz = linspace(18, 44, 401)';
    f    = fGHz * 1e9;
    k0   = 2 * pi * f / c0;

    % =================================================================
    % PART 1: TMM Forward Model Validation
    % =================================================================
    fprintf('=== Part 1: TMM Forward Model ===\n\n');

    % Typical wood parameters at microwave frequencies
    % Earlywood: less dense, more moisture -> higher eps', higher loss
    % Latewood:  denser, less moisture   -> lower eps', lower loss

    % Reference: from literature, dry wood eps' ~ 2-3, water eps' ~ 40-80
    % Moist wood effective: eps' ~ 3-8, tan(delta) ~ 0.05-0.20

    earlywood = struct('eps_p', 5.0, 'eps_pp', 1.0, 'd_mm', 3.0);
    latewood  = struct('eps_p', 3.0, 'eps_pp', 0.5, 'd_mm', 2.0);

    d_total = 0.120;  % 120 mm total thickness

    % Compute number of layer pairs
    d_pair_mm = earlywood.d_mm + latewood.d_mm;
    N_pairs = round(d_total * 1000 / d_pair_mm);
    N_layers = 2 * N_pairs;
    fprintf('Total thickness: %.0f mm\n', d_total*1000);
    fprintf('Layer pair (E+L): %.1f mm, N_pairs = %d, N_layers = %d\n', ...
        d_pair_mm, N_pairs, N_layers);
    fprintf('Earlywood: eps=%.2f-j%.2f, d=%.1fmm\n', ...
        earlywood.eps_p, earlywood.eps_pp, earlywood.d_mm);
    fprintf('Latewood:  eps=%.2f-j%.2f, d=%.1fmm\n', ...
        latewood.eps_p, latewood.eps_pp, latewood.d_mm);

    % Compute TMM S-parameters
    [S11_tmm, S21_tmm] = tmm_layered(earlywood, latewood, N_pairs, f, c0, Z0);

    % Compute equivalent uniform medium for Fabry-Perot comparison
    % Volume-weighted average eps*
    frac_E = earlywood.d_mm / d_pair_mm;
    frac_L = latewood.d_mm / d_pair_mm;
    eps_avg = (frac_E * (earlywood.eps_p - 1j*earlywood.eps_pp) + ...
               frac_L * (latewood.eps_p - 1j*latewood.eps_pp));
    eps_avg_p  = real(eps_avg);
    eps_avg_pp = -imag(eps_avg);

    % Fabry-Perot with averaged medium
    S21_fp = fabry_perot_uniform(eps_avg_p, eps_avg_pp, d_total, f, c0);

    % ---- Plot comparison ----
    fig1 = figure('Color', 'w', 'Position', [50, 100, 1100, 700]);

    subplot(2,2,1);
    plot(fGHz, 20*log10(abs(S21_tmm)), 'b-', 'LineWidth', 1.5); hold on;
    plot(fGHz, 20*log10(abs(S21_fp)), 'r--', 'LineWidth', 1.5);
    xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 12);
    ylabel('|S_{21}| (dB)', 'FontName', fontName, 'FontSize', 12);
    legend({'TMM (layered)', 'Fabry-Perot (uniform)'}, ...
        'FontName', fontName, 'FontSize', 10, 'Location', 'southwest');
    title('|S_{21}|: Layered vs Uniform', ...
        'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', 11, 'Box', 'on'); grid on;

    subplot(2,2,2);
    res = 20*log10(abs(S21_tmm)) - 20*log10(abs(S21_fp));
    plot(fGHz, res, 'k-', 'LineWidth', 1.5); hold on;
    yline(0, 'k--');
    xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 12);
    ylabel('Residual (dB)', 'FontName', fontName, 'FontSize', 12);
    title(sprintf('TMM - FP Residual (RMS=%.4f dB)', my_rms(res)), ...
        'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', 11, 'Box', 'on'); grid on;

    subplot(2,2,3);
    eps_eff_tmm = -1 ./ (1j * k0 * d_total) .* log(S21_tmm);
    plot(fGHz, real(eps_eff_tmm), 'b-', 'LineWidth', 1.5); hold on;
    plot(fGHz, [-imag(eps_eff_tmm)], 'r-', 'LineWidth', 1.5);
    yline(eps_avg_p, 'b--'); yline(eps_avg_pp, 'r--');
    xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 12);
    ylabel('Effective \epsilon^*', 'FontName', fontName, 'FontSize', 12);
    legend({'\epsilon''', '\epsilon''''', 'Avg \epsilon''', 'Avg \epsilon'''''}, ...
        'FontName', fontName, 'FontSize', 9, 'Location', 'best');
    title('Effective Permittivity from S_{21}', ...
        'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', 11, 'Box', 'on'); grid on;

    subplot(2,2,4);
    plot(fGHz, 20*log10(abs(S11_tmm)), 'b-', 'LineWidth', 1.5);
    xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 12);
    ylabel('|S_{11}| (dB)', 'FontName', fontName, 'FontSize', 12);
    title('|S_{11}|: Reflection (TMM)', ...
        'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', 11, 'Box', 'on'); grid on;

    sgtitle(sprintf(['Layered Wood TMM: %d layers, E(%.1fmm)/L(%.1fmm), ' ...
        '\\Delta\\epsilon''=%.1f'], ...
        N_layers, earlywood.d_mm, latewood.d_mm, ...
        earlywood.eps_p - latewood.eps_p), ...
        'FontName', fontName, 'FontSize', 14, 'FontWeight', 'bold');
    exportgraphics(fig1, fullfile(outDir, 'Exp6_TMM_vs_FP.png'), 'Resolution', 300);
    fprintf('  Saved: Exp6_TMM_vs_FP.png\n');

    % =================================================================
    % PART 2: Dielectric Contrast Sweep
    % =================================================================
    fprintf('\n=== Part 2: Dielectric Contrast Sensitivity ===\n\n');

    % Sweep earlywood-latewood contrast: how different must they be
    % for the layered structure to matter?
    contrast_vals = 0.5:0.5:4.0;  % delta_eps' between E and L
    my_rms_tmm_fp = zeros(size(contrast_vals));

    fprintf('%-15s %10s\n', 'Contrast', 'RMS(dB)');
    fprintf('%s\n', repmat('-', 1, 30));

    for iC = 1:length(contrast_vals)
        dc = contrast_vals(iC);
        ew = earlywood;
        lw = latewood;
        % Vary contrast: keep average constant, split around it
        avg_p = 4.0;
        ew.eps_p = avg_p + dc/2;
        lw.eps_p = max(1.5, avg_p - dc/2);

        [~, S21_t] = tmm_layered(ew, lw, N_pairs, f, c0, Z0);
        eps_avg2 = (frac_E*(ew.eps_p - 1j*ew.eps_pp) + ...
                    frac_L*(lw.eps_p - 1j*lw.eps_pp));
        S21_f = fabry_perot_uniform(real(eps_avg2), -imag(eps_avg2), ...
                                    d_total, f, c0);
        my_rms_tmm_fp(iC) = my_rms(20*log10(abs(S21_t)) - 20*log10(abs(S21_f)));
        fprintf('%-15.1f %10.3f\n', dc, my_rms_tmm_fp(iC));
    end

    fig2 = figure('Color', 'w', 'Position', [50, 100, 900, 400]);

    subplot(1,2,1);
    plot(contrast_vals, my_rms_tmm_fp, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
    xlabel('\Delta\epsilon'' (Earlywood - Latewood)', ...
        'FontName', fontName, 'FontSize', 12);
    ylabel('RMS Residual (dB)', 'FontName', fontName, 'FontSize', 12);
    title('Layered vs Uniform: Effect of Contrast', ...
        'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', 11, 'Box', 'on'); grid on;
    yline(0.1, 'r--', 'LineWidth', 1.0);
    text(3.0, 0.12, 'Detection threshold ~0.1 dB', ...
        'FontName', fontName, 'FontSize', 9, 'Color', 'r');

    % Show S21 for different contrasts
    subplot(1,2,2);
    cmap = jet(length(contrast_vals));
    hold on;
    for iC = 1:2:length(contrast_vals)
        dc = contrast_vals(iC);
        ew = earlywood; lw = latewood;
        avg_p = 4.0;
        ew.eps_p = avg_p + dc/2;
        lw.eps_p = max(1.5, avg_p - dc/2);
        [~, S21_t] = tmm_layered(ew, lw, N_pairs, f, c0, Z0);
        plot(fGHz, 20*log10(abs(S21_t)), 'Color', cmap(iC,:), 'LineWidth', 1.2);
    end
    xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 12);
    ylabel('|S_{21}| (dB)', 'FontName', fontName, 'FontSize', 12);
    title('|S_{21}| at Different Contrasts', ...
        'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', 11, 'Box', 'on'); grid on;
    cb = colorbar; colormap(jet);
    cb.Label.String = '\Delta\epsilon''';
    cb.Label.FontName = fontName;

    sgtitle('Layered Structure: When Does It Matter?', ...
        'FontName', fontName, 'FontSize', 14, 'FontWeight', 'bold');
    exportgraphics(fig2, fullfile(outDir, 'Exp6_Contrast_Sweep.png'), 'Resolution', 300);
    fprintf('  Saved: Exp6_Contrast_Sweep.png\n');

    % =================================================================
    % PART 3: Parameter Identifiability (Reduced Model)
    % =================================================================
    fprintf('\n=== Part 3: Reduced-Model Identifiability ===\n\n');

    % We cannot identify 6 params. Use physically-motivated reductions:
    %
    % Assumption 1: Same loss tangent for E and L (both wood, similar moisture)
    %   -> tan_delta_E = tan_delta_L = tan_delta
    %
    % Assumption 2: Thickness ratio known from annual ring measurement
    %   -> d_E / d_L fixed (e.g., 3:2 typical for temperate wood)
    %
    % Then: X = [eps_p_E, eps_p_L, tan_delta, d_pair]
    % 4 params — same as uniform model, but with layered physics.

    % True parameters for reduced model
    truth_red = [5.0, 3.0, 0.15, 5.0];  % [eps_E, eps_L, tand, d_pair_mm]
    d_E_true = 3.0;   % mm
    d_L_true = 2.0;   % mm

    % Generate TMM data
    ew_true = struct('eps_p', truth_red(1), 'eps_pp', truth_red(1)*truth_red(3), ...
                     'd_mm', d_E_true);
    lw_true = struct('eps_p', truth_red(2), 'eps_pp', truth_red(2)*truth_red(3), ...
                     'd_mm', d_L_true);
    N_pairs_true = round(d_total * 1000 / (d_E_true + d_L_true));

    [~, S21_true] = tmm_layered(ew_true, lw_true, N_pairs_true, f, c0, Z0);
    mag_true_dB = 20*log10(abs(S21_true));

    % 1D sweeps
    epsE_vals  = linspace(3.0, 7.0, 60);
    epsL_vals  = linspace(2.0, 5.0, 60);
    tand_vals  = linspace(0.05, 0.30, 60);
    dpair_vals = linspace(3.0, 7.0, 60);

    cE   = sweep_tmm_1d('eps_E', epsE_vals, truth_red, d_E_true, d_L_true, ...
                         N_pairs_true, d_total, f, c0, Z0, mag_true_dB);
    cL   = sweep_tmm_1d('eps_L', epsL_vals, truth_red, d_E_true, d_L_true, ...
                         N_pairs_true, d_total, f, c0, Z0, mag_true_dB);
    cT   = sweep_tmm_1d('tand',  tand_vals, truth_red, d_E_true, d_L_true, ...
                         N_pairs_true, d_total, f, c0, Z0, mag_true_dB);
    cD   = sweep_tmm_1d('d_pair', dpair_vals, truth_red, d_E_true, d_L_true, ...
                         N_pairs_true, d_total, f, c0, Z0, mag_true_dB);

    sens = zeros(4,1);
    all_c = {cE, cL, cT, cD};
    for k = 1:4, sens(k) = (max(all_c{k}) - min(all_c{k})) / mean(all_c{k}); end

    fig3 = figure('Color', 'w', 'Position', [50, 100, 1000, 700]);

    param_info = { ...
        {epsE_vals, cE, '\epsilon''_E (earlywood)', truth_red(1)}, ...
        {epsL_vals, cL, '\epsilon''_L (latewood)',  truth_red(2)}, ...
        {tand_vals, cT, 'tan\delta (shared)',       truth_red(3)}, ...
        {dpair_vals, cD, 'd_{pair} (E+L, mm)',      truth_red(4)}  ...
    };

    for k = 1:4
        subplot(2, 2, k);
        plot(param_info{k}{1}, param_info{k}{2}, 'b-', 'LineWidth', 1.5); hold on;
        xline(param_info{k}{4}, 'r--', 'LineWidth', 1.2);
        xlabel(param_info{k}{3}, 'FontName', fontName, 'FontSize', 12);
        ylabel('RMSE (dB)', 'FontName', fontName, 'FontSize', 12);
        title(sprintf('Cost vs %s  (sens=%.2f)', param_info{k}{3}, sens(k)), ...
            'FontName', fontName, 'FontSize', 11, 'FontWeight', 'bold');
        set(gca, 'FontName', fontName, 'FontSize', 10, 'Box', 'on'); grid on;
    end

    sgtitle('TMM Reduced Model: 4-Parameter Sensitivity', ...
        'FontName', fontName, 'FontSize', 14, 'FontWeight', 'bold');
    exportgraphics(fig3, fullfile(outDir, 'Exp6_TMM_Sensitivity.png'), 'Resolution', 300);
    fprintf('  Saved: Exp6_TMM_Sensitivity.png\n');

    % =================================================================
    % PART 4: Recovery Test with Reduced TMM Model
    % =================================================================
    fprintf('\n=== Part 4: TMM Recovery Test ===\n\n');

    lb_red = [2.0, 1.5, 0.02, 3.0];
    ub_red = [7.0, 5.0, 0.40, 8.0];
    nTrials = 8;
    noise_dB = 0.1;

    fprintf('%-5s %10s %10s %10s %10s\n', 'Trial', 'd(eps_E)', 'd(eps_L)', 'd(tand)', 'd(dpair)');
    fprintf('%s\n', repmat('-', 1, 55));

    rec = zeros(nTrials, 4);
    for iT = 1:nTrials
        rng(200 + iT);
        mag_n = mag_true_dB + noise_dB * randn(size(mag_true_dB));
        T_noisy = 10.^(mag_n/20) .* exp(1j * angle(S21_true));

        [prec, ~] = tmm_invert_red(T_noisy, f, lb_red, ub_red, d_total, ...
                                   40, 300, c0, Z0);
        err = prec - truth_red;
        rec(iT, :) = prec;
        fprintf('%5d %+10.3f %+10.3f %+10.3f %+10.2f\n', iT, err(1), err(2), err(3), err(4));
    end
    err_mean = mean(rec, 1) - truth_red;
    err_std  = std(rec, 0, 1);
    fprintf('Mean:  %+10.3f %+10.3f %+10.3f %+10.2f\n', ...
        err_mean(1), err_mean(2), err_mean(3), err_mean(4));
    fprintf('Std:   %10.3f %10.3f %10.3f %10.2f\n', ...
        err_std(1), err_std(2), err_std(3), err_std(4));

    % ---- Plot recovery ----
    rng(42);
    mag_n = mag_true_dB + noise_dB * randn(size(mag_true_dB));
    T_noisy = 10.^(mag_n/20) .* exp(1j * angle(S21_true));
    [p_best, loss_best] = tmm_invert_red(T_noisy, f, lb_red, ub_red, d_total, ...
                                         50, 400, c0, Z0);

    ew_best = struct('eps_p', p_best(1), 'eps_pp', p_best(1)*p_best(3), ...
                     'd_mm', d_E_true);
    lw_best = struct('eps_p', p_best(2), 'eps_pp', p_best(2)*p_best(3), ...
                     'd_mm', d_L_true);
    N_best = round(d_total * 1000 / (d_E_true + d_L_true));
    [~, S21_best] = tmm_layered(ew_best, lw_best, N_best, f, c0, Z0);

    % Also fit with uniform FP for comparison
    [pFP, ~] = debye3_invert(T_noisy, f, [1.5 0.5 1e9], [5 6 80e9], ...
                              d_total, 50, 400, c0);
    T_FP = debye_T(pFP(1), pFP(2), pFP(3), d_total, f, c0);

    fig4 = figure('Color', 'w', 'Position', [50, 100, 1100, 450]);

    subplot(1,2,1);
    plot(fGHz, mag_n, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 3); hold on;
    plot(fGHz, 20*log10(abs(S21_true)), 'k-', 'LineWidth', 1.8);
    plot(fGHz, 20*log10(abs(S21_best)), 'b--', 'LineWidth', 1.5);
    plot(fGHz, 20*log10(abs(T_FP)), 'r-.', 'LineWidth', 1.5);
    xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 12);
    ylabel('|S_{21}| (dB)', 'FontName', fontName, 'FontSize', 12);
    legend({'Noisy data', 'TMM truth', 'TMM fit', 'Uniform FP fit'}, ...
        'FontName', fontName, 'FontSize', 9, 'Location', 'southwest');
    title('|S_{21}| Fits: TMM vs Uniform FP', ...
        'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', 11, 'Box', 'on'); grid on;

    subplot(1,2,2);
    res_tmm = 20*log10(abs(S21_best)) - 20*log10(abs(S21_true));
    res_fp  = 20*log10(abs(T_FP)) - 20*log10(abs(S21_true));
    plot(fGHz, res_tmm, 'b-', 'LineWidth', 1.5); hold on;
    plot(fGHz, res_fp, 'r-', 'LineWidth', 1.5);
    yline(0, 'k--');
    xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 12);
    ylabel('Residual (dB)', 'FontName', fontName, 'FontSize', 12);
    legend({sprintf('TMM (RMS=%.4f dB)', my_rms(res_tmm)), ...
            sprintf('FP (RMS=%.4f dB)', my_rms(res_fp))}, ...
        'FontName', fontName, 'FontSize', 10, 'Location', 'best');
    title('Fitting Residual Comparison', ...
        'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', 11, 'Box', 'on'); grid on;

    sgtitle('TMM Layered Model Recovery: 4-Parameter Reduced Form', ...
        'FontName', fontName, 'FontSize', 14, 'FontWeight', 'bold');
    exportgraphics(fig4, fullfile(outDir, 'Exp6_TMM_Recovery.png'), 'Resolution', 300);
    fprintf('\n  Saved: Exp6_TMM_Recovery.png\n');

    % ---- Summary ----
    fprintf('\n=== TMM Key Findings ===\n');
    fprintf('Contrast threshold for noticeable layering: ~%.1f in Delta_eps\n', ...
        contrast_vals(find(my_rms_tmm_fp > 0.1, 1, 'first')));
    fprintf('TMM RMS: %.4f dB, FP RMS: %.4f dB\n', my_rms(res_tmm), my_rms(res_fp));
    fprintf('\n=== Exp6 complete. Results in: %s ===\n', outDir);
end

%% ==================== TMM FORWARD MODEL ====================

function [S11, S21] = tmm_layered(earlywood, latewood, N_pairs, f, c0, Z0)
    % Build alternating E-L-E-L-...-E-L stack using Chebyshev method.
    % Vectorized over frequency for speed (no per-frequency loop).

    eps_E = earlywood.eps_p - 1j * earlywood.eps_pp;
    eps_L = latewood.eps_p  - 1j * latewood.eps_pp;
    d_E   = earlywood.d_mm * 1e-3;
    d_L   = latewood.d_mm  * 1e-3;

    n_E = sqrt(eps_E); n_L = sqrt(eps_L);
    Z_E = Z0 / n_E;  Z_L = Z0 / n_L;

    k0 = (2 * pi * f / c0);  % N×1

    % Phase delays per layer
    phi_E = k0 * n_E * d_E;
    phi_L = k0 * n_L * d_L;

    % M_pair = M_E * M_L for each frequency
    % M_E = [cos(phiE), j*ZE*sin(phiE); j/ZE*sin(phiE), cos(phiE)]
    % M_L = [cos(phiL), j*ZL*sin(phiL); j/ZL*sin(phiL), cos(phiL)]

    cE = cos(phi_E); sE = sin(phi_E);
    cL = cos(phi_L); sL = sin(phi_L);

    % M_pair = [A B; C D] (element-wise for all frequencies)
    A = cE.*cL - (Z_E/Z_L).*sE.*sL;
    B = 1j*(Z_E.*sE.*cL + Z_L.*cE.*sL);
    C = 1j*(1./Z_E.*sE.*cL + 1./Z_L.*cE.*sL);
    D = cE.*cL - (Z_L/Z_E).*sE.*sL;

    % Use Chebyshev identity: M_pair^N = U_{N-1} * M_pair - U_{N-2} * I
    % where U_n(cos(theta)) = sin((n+1)*theta)/sin(theta)
    s = (A + D) / 2;  % = trace(M)/2

    % For |s| <= 1 (propagating): theta = acos(s)
    % For |s| > 1  (evanescent): x = acosh(s)
    mask_prop = abs(s) <= 1;

    % Compute Chebyshev polynomials (guard against NaN when s -> +/-1)
    theta = acos(s);
    sin_th = sin(theta);
    near_degen = abs(sin_th) < 1e-7;

    U_Nm1 = zeros(size(s));
    U_Nm2 = zeros(size(s));

    if any(near_degen)
        sgn = sign(real(s(near_degen)));
        U_Nm1(near_degen) = N_pairs * sgn.^(N_pairs-1);
        U_Nm2(near_degen) = (N_pairs-1) * sgn.^(N_pairs-2);
    end
    ok = ~near_degen;
    U_Nm1(ok) = sin(N_pairs * theta(ok)) ./ sin_th(ok);
    U_Nm2(ok) = sin((N_pairs-1) * theta(ok)) ./ sin_th(ok);

    % M_total = U_Nm1 * M_pair - U_Nm2 * I
    At = U_Nm1 .* A - U_Nm2;
    Bt = U_Nm1 .* B;
    Ct = U_Nm1 .* C;
    Dt = U_Nm1 .* D - U_Nm2;

    % Convert ABCD to S-parameters (both sides Z0)
    denom = At + Bt/Z0 + Ct*Z0 + Dt;
    S11 = (At + Bt/Z0 - Ct*Z0 - Dt) ./ denom;
    S21 = 2 ./ denom;
end

function M = layer_matrix(eps_c, d, k0, Z0)
    n_complex = sqrt(eps_c);
    phi = k0 * n_complex * d;
    Z = Z0 / n_complex;

    M = [cos(phi), 1j*Z*sin(phi);
         1j/Z*sin(phi), cos(phi)];
end

%% ==================== UNIFORM FABRY-PEROT ====================

function T = fabry_perot_uniform(eps_p, eps_pp, d, f, c0)
    eps_c = eps_p - 1j * eps_pp;
    gamma = 1j * (2*pi*f) / c0 .* sqrt(eps_c);
    R = (sqrt(eps_c) - 1) / (sqrt(eps_c) + 1);
    E2 = exp(-2 * gamma * d);
    T = exp(-gamma * d) .* (1 - R^2) ./ (1 - R^2 * E2);
end

%% ==================== TMM 1D SWEEP (REDUCED MODEL) ====================

function cost = sweep_tmm_1d(name, vals, truth_red, d_E_mm, d_L_mm, ...
        N_pairs, d_total, f, c0, Z0, mag_true_dB)
    Nv = length(vals);
    cost = zeros(Nv, 1);
    for i = 1:Nv
        p = truth_red;
        switch name
            case 'eps_E',  p(1) = vals(i);
            case 'eps_L',  p(2) = vals(i);
            case 'tand',   p(3) = vals(i);
            case 'd_pair', p(4) = vals(i);
        end
        ew = struct('eps_p', p(1), 'eps_pp', p(1)*p(3), 'd_mm', d_E_mm);
        lw = struct('eps_p', p(2), 'eps_pp', p(2)*p(3), 'd_mm', d_L_mm);
        Np = round(d_total * 1000 / (d_E_mm + d_L_mm));
        [~, S21] = tmm_layered(ew, lw, Np, f, c0, Z0);
        cost(i) = sqrt(mean((20*log10(abs(S21)) - mag_true_dB).^2));
    end
end

%% ==================== TMM INVERSION (REDUCED 4-PARAM) ====================

function [best_x, loss_curve] = tmm_invert_red(T_meas, f, lb, ub, d_total, ...
        popSize, maxIter, c0, Z0)
    D = 4;
    mag_meas_dB = 20 * log10(abs(T_meas));

    X = repmat(lb, popSize, 1) + lhsdesign(popSize, D) .* repmat(ub - lb, popSize, 1);

    cost = zeros(popSize, 1);
    for i = 1:popSize
        cost(i) = cost_tmm_red(X(i,:), f, mag_meas_dB, d_total, c0, Z0);
    end
    [best_cost, idx] = min(cost);
    best_x = X(idx, :);
    loss_curve = zeros(maxIter, 1);

    F = 0.5; CR = 0.7;
    for iter = 1:maxIter
        for i = 1:popSize
            candidates = setdiff(1:popSize, i);
            r = candidates(randperm(length(candidates), 3));
            v = X(r(1), :) + F * (X(r(2), :) - X(r(3), :));
            j_rand = randi(D);
            u = X(i, :);
            for j = 1:D
                if rand < CR || j == j_rand, u(j) = v(j); end
            end
            u = max(lb, min(ub, u));
            cost_u = cost_tmm_red(u, f, mag_meas_dB, d_total, c0, Z0);
            if cost_u < cost(i)
                X(i, :) = u;
                cost(i) = cost_u;
                if cost_u < best_cost
                    best_cost = cost_u;
                    best_x = u;
                end
            end
        end
        loss_curve(iter) = best_cost;
    end
end

function J = cost_tmm_red(x, f, mag_meas_dB, d_total, c0, Z0)
    % x = [eps_p_E, eps_p_L, tand_shared, d_pair_mm]
    % Fixed: thickness ratio d_E:d_L = 3:2
    d_E_mm = x(4) * 0.6;
    d_L_mm = x(4) * 0.4;
    ew = struct('eps_p', x(1), 'eps_pp', x(1)*x(3), 'd_mm', d_E_mm);
    lw = struct('eps_p', x(2), 'eps_pp', x(2)*x(3), 'd_mm', d_L_mm);
    Np = round(d_total * 1000 / (d_E_mm + d_L_mm));
    [~, S21] = tmm_layered(ew, lw, Np, f, c0, Z0);
    mag_err = 20*log10(abs(S21)) - mag_meas_dB;
    Lm = mag_err.^2;
    mask = abs(mag_err) > 1.0;
    Lm(mask) = 2.0 * abs(mag_err(mask)) - 1.0;
    J = mean(Lm);
end

%% ==================== DEBYE (for comparison) ====================

function T = debye_T(eps_inf, delta_eps, f_relax, d, f, c0)
    eps_c = eps_inf + delta_eps ./ (1 + 1j * f / f_relax);
    gamma = 1j * (2*pi*f) / c0 .* sqrt(eps_c);
    R = (sqrt(eps_c) - 1) ./ (sqrt(eps_c) + 1);
    E2 = exp(-2 * gamma * d);
    T = exp(-gamma * d) .* (1 - R.^2) ./ (1 - R.^2 .* E2);
end

function [best_x, loss_curve] = debye3_invert(T_meas, f, lb, ub, d_fixed, ...
        popSize, maxIter, c0)
    D = 3;
    mag_meas_dB = 20 * log10(abs(T_meas));
    X = repmat(lb, popSize, 1) + lhsdesign(popSize, D) .* repmat(ub - lb, popSize, 1);
    cost = zeros(popSize, 1);
    for i = 1:popSize
        T = debye_T(X(i,1), X(i,2), X(i,3), d_fixed, f, c0);
        mag_err = 20*log10(abs(T)) - mag_meas_dB;
        Lm = mag_err.^2; mask = abs(mag_err) > 1.0;
        Lm(mask) = 2.0 * abs(mag_err(mask)) - 1.0;
        cost(i) = mean(Lm);
    end
    [best_cost, idx] = min(cost);
    best_x = X(idx, :);
    loss_curve = zeros(maxIter, 1);
    F = 0.5; CR = 0.7;
    for iter = 1:maxIter
        for i = 1:popSize
            candidates = setdiff(1:popSize, i);
            r = candidates(randperm(length(candidates), 3));
            v = X(r(1), :) + F * (X(r(2), :) - X(r(3), :));
            j_rand = randi(D);
            u = X(i, :);
            for j = 1:D
                if rand < CR || j == j_rand, u(j) = v(j); end
            end
            u = max(lb, min(ub, u));
            T = debye_T(u(1), u(2), u(3), d_fixed, f, c0);
            mag_err = 20*log10(abs(T)) - mag_meas_dB;
            Lm = mag_err.^2; mask = abs(mag_err) > 1.0;
            Lm(mask) = 2.0 * abs(mag_err(mask)) - 1.0;
            cost_u = mean(Lm);
            if cost_u < cost(i)
                X(i, :) = u; cost(i) = cost_u;
                if cost_u < best_cost
                    best_cost = cost_u; best_x = u;
                end
            end
        end
        loss_curve(iter) = best_cost;
    end
end

function r = my_rms(x)
    r = sqrt(mean(x(:).^2, 'omitnan'));
end
