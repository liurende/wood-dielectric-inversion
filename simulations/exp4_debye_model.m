%% Experiment 4: Debye Relaxation Model for Wood
% Validate the Debye model identifiability and compare with power-law.
% Key question: can we recover Debye parameters from S21 data?

function exp4_debye_model()
    clc;

    outDir = fullfile(pwd, 'simulations', 'results');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    fontName = 'Times New Roman';
    c0 = 299792458;
    fGHz = linspace(18, 44, 201)';
    f    = fGHz * 1e9;

    % ---- Reference Debye parameters (moist wood) ----
    p_true = struct('eps_inf', 2.5, 'delta_eps', 2.0, ...
                    'f_relax', 15e9, 'd', 0.120);
    fprintf('Reference: eps_inf=%.1f, delta_eps=%.1f, f_relax=%.0fGHz, d=%.0fmm\n', ...
        p_true.eps_inf, p_true.delta_eps, p_true.f_relax/1e9, p_true.d*1000);

    % ---- Compute effective permittivity at 30 GHz (for comparison) ----
    eps_at_ref = debye_eps(p_true.eps_inf, p_true.delta_eps, p_true.f_relax, 30e9);
    fprintf('Effective at 30 GHz: epsilon* = %.3f - j%.4f  (tan_delta=%.4f)\n\n', ...
        real(eps_at_ref), -imag(eps_at_ref), -imag(eps_at_ref)/real(eps_at_ref));

    % =================================================================
    % PART 1: 1D Parameter Sensitivity
    % =================================================================
    fprintf('Part 1: 1D sensitivity sweeps...\n');

    eps_inf_vals  = linspace(1.5, 5.0, 80);
    delta_eps_vals = linspace(0.5, 6.0, 80);
    frelax_vals   = logspace(log10(2e9), log10(80e9), 80);
    d_vals        = linspace(0.05, 0.40, 80);

    T_true = debye_T(p_true.eps_inf, p_true.delta_eps, p_true.f_relax, ...
                     p_true.d, f, c0);
    mag_true_dB = 20 * log10(abs(T_true));

    c1 = sweep_debye_1d('eps_inf',  eps_inf_vals,  p_true, f, c0, mag_true_dB);
    c2 = sweep_debye_1d('delta_eps', delta_eps_vals, p_true, f, c0, mag_true_dB);
    c3 = sweep_debye_1d('f_relax',  frelax_vals,   p_true, f, c0, mag_true_dB);
    c4 = sweep_debye_1d('d',        d_vals,        p_true, f, c0, mag_true_dB);

    % Normalize to min-max for fair comparison
    c_all = {c1, c2, c3, c4};
    curv = zeros(4,1);
    for k = 1:4, curv(k) = (max(c_all{k}) - min(c_all{k})) / mean(c_all{k}); end

    fig1 = figure('Color', 'w', 'Position', [50, 100, 1000, 700]);

    subplot(2,2,1);
    plot(eps_inf_vals, c1, 'b-', 'LineWidth', 1.5); hold on;
    xline(p_true.eps_inf, 'r--', 'LineWidth', 1.2);
    xlabel('\epsilon_\infty', 'FontName', fontName, 'FontSize', 14);
    ylabel('RMSE (dB)', 'FontName', fontName, 'FontSize', 14);
    title(sprintf('Cost vs \\epsilon_\\infty  (sens=%.2f)', curv(1)), ...
        'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', 12, 'Box', 'on'); grid on;

    subplot(2,2,2);
    plot(delta_eps_vals, c2, 'b-', 'LineWidth', 1.5); hold on;
    xline(p_true.delta_eps, 'r--', 'LineWidth', 1.2);
    xlabel('\Delta\epsilon', 'FontName', fontName, 'FontSize', 14);
    ylabel('RMSE (dB)', 'FontName', fontName, 'FontSize', 14);
    title(sprintf('Cost vs \\Delta\\epsilon  (sens=%.2f)', curv(2)), ...
        'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', 12, 'Box', 'on'); grid on;

    subplot(2,2,3);
    semilogx(frelax_vals/1e9, c3, 'b-', 'LineWidth', 1.5); hold on;
    xline(p_true.f_relax/1e9, 'r--', 'LineWidth', 1.2);
    xlabel('f_{relax} (GHz)', 'FontName', fontName, 'FontSize', 14);
    ylabel('RMSE (dB)', 'FontName', fontName, 'FontSize', 14);
    title(sprintf('Cost vs f_{relax}  (sens=%.2f)', curv(3)), ...
        'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', 12, 'Box', 'on'); grid on;

    subplot(2,2,4);
    plot(d_vals*1000, c4, 'b-', 'LineWidth', 1.5); hold on;
    xline(p_true.d*1000, 'r--', 'LineWidth', 1.2);
    xlabel('d (mm)', 'FontName', fontName, 'FontSize', 14);
    ylabel('RMSE (dB)', 'FontName', fontName, 'FontSize', 14);
    title(sprintf('Cost vs d  (sens=%.2f)', curv(4)), ...
        'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', 12, 'Box', 'on'); grid on;

    sgtitle('Debye Model: 1D Parameter Sensitivity', ...
        'FontName', fontName, 'FontSize', 16, 'FontWeight', 'bold');
    exportgraphics(fig1, fullfile(outDir, 'Exp4_Debye_1D_Sensitivity.png'), 'Resolution', 300);
    fprintf('  Saved: Exp4_Debye_1D_Sensitivity.png\n');

    % =================================================================
    % PART 2: 2D Coupling Analysis
    % =================================================================
    fprintf('Part 2: 2D coupling contours...\n');

    np = 50;  % coarser grid for speed
    ei_vals = linspace(1.5, 5.0, np);
    de_vals = linspace(0.5, 6.0, np);
    fr_vals = logspace(log10(2e9), log10(80e9), np);
    dd_vals = linspace(0.05, 0.40, np);

    [C_inf_de,  ~, ~] = sweep_debye_2d('eps_inf', ei_vals, 'delta_eps', de_vals, ...
        p_true, f, c0, mag_true_dB);
    [C_inf_fr,  ~, ~] = sweep_debye_2d('eps_inf', ei_vals, 'f_relax', fr_vals, ...
        p_true, f, c0, mag_true_dB);
    [C_de_fr,   ~, ~] = sweep_debye_2d('delta_eps', de_vals, 'f_relax', fr_vals, ...
        p_true, f, c0, mag_true_dB);
    [C_inf_d,   ~, ~] = sweep_debye_2d('eps_inf', ei_vals, 'd', dd_vals, ...
        p_true, f, c0, mag_true_dB);

    ar_inf_de = coupling_ratio(C_inf_de, ei_vals, de_vals);
    ar_inf_fr = coupling_ratio(C_inf_fr, ei_vals, fr_vals);
    ar_de_fr  = coupling_ratio(C_de_fr,  de_vals, fr_vals);
    ar_inf_d  = coupling_ratio(C_inf_d,  ei_vals, dd_vals);

    fig2 = figure('Color', 'w', 'Position', [50, 100, 1100, 850]);

    pairs = { ...
        {ei_vals, de_vals, C_inf_de, '\epsilon_\infty', '\Delta\epsilon', ar_inf_de}, ...
        {ei_vals, fr_vals, C_inf_fr, '\epsilon_\infty', 'f_{relax} (GHz)', ar_inf_fr}, ...
        {de_vals, fr_vals, C_de_fr,  '\Delta\epsilon', 'f_{relax} (GHz)', ar_de_fr}, ...
        {ei_vals, dd_vals, C_inf_d,  '\epsilon_\infty', 'd (mm)', ar_inf_d}  ...
    };

    for k = 1:4
        subplot(2, 2, k);
        xv = pairs{k}{1}; yv = pairs{k}{2}; C = pairs{k}{3};
        if k == 2 || k == 3
            xv_disp = log10(xv);
            xtick_labels = {'2','5','10','20','50'};
        else
            xv_disp = xv;
        end
        if k == 2 || k == 3
            yv_disp = log10(yv);
        else
            yv_disp = yv;
        end
        if k == 4, yv_disp = yv * 1000; end

        contourf(xv_disp, yv_disp, C', 30, 'LineStyle', 'none');
        hold on;
        % Mark true value
        switch k
            case 1
                plot(p_true.eps_inf, p_true.delta_eps, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
            case 2
                plot(log10(p_true.f_relax), log10(p_true.eps_inf), 'rx', ...
                    'MarkerSize', 12, 'LineWidth', 2);
            case 3
                plot(log10(p_true.f_relax), log10(p_true.delta_eps), 'rx', ...
                    'MarkerSize', 12, 'LineWidth', 2);
            case 4
                plot(p_true.eps_inf, p_true.d*1000, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
        end
        xlabel(pairs{k}{4}, 'FontName', fontName, 'FontSize', 13);
        ylabel(pairs{k}{5}, 'FontName', fontName, 'FontSize', 13);
        if k == 2 || k == 3
            set(gca, 'XTick', log10([2,5,10,20,50]), 'XTickLabel', {'2','5','10','20','50'});
            set(gca, 'YTick', log10([2,5,10,20,50]));
            set(gca, 'YTickLabel', {});
        end
        title(sprintf('%s vs %s  (AR=%.2f)', pairs{k}{4}, pairs{k}{5}, pairs{k}{6}), ...
            'FontName', fontName, 'FontSize', 11, 'FontWeight', 'bold');
        colorbar; colormap(jet);
        set(gca, 'FontName', fontName, 'FontSize', 10, 'Box', 'on');
    end

    sgtitle('Debye Model: 2D Parameter Coupling', ...
        'FontName', fontName, 'FontSize', 16, 'FontWeight', 'bold');
    exportgraphics(fig2, fullfile(outDir, 'Exp4_Debye_2D_Coupling.png'), 'Resolution', 300);
    fprintf('  Saved: Exp4_Debye_2D_Coupling.png\n');

    % =================================================================
    % PART 3: Recovery test (Debye inversion on Debye data)
    % =================================================================
    fprintf('Part 3: Debye recovery test...\n');

    % Test 3 scenarios: f_relax inside band, below band, above band
    scenarios = [
        2.5, 2.0, 15.0, 0.120;   % f_relax inside band (18-44 GHz)
        2.5, 2.0,  6.0, 0.120;   % f_relax below band
        2.5, 2.0, 60.0, 0.120;   % f_relax above band
    ];
    scNames = {'f_{relax} in band (15 GHz)', ...
               'f_{relax} below band (6 GHz)', ...
               'f_{relax} above band (60 GHz)'};

    lb_debye = [1.5, 0.5,  1e9, 0.05];
    ub_debye = [5.0, 6.0, 100e9, 0.40];
    popSize = 50;
    maxIter = 300;
    nTrials = 8;
    noise_mag_dB = 0.1;
    noise_phase_deg = 1.0;

    fprintf('%-25s %10s %10s %10s %10s\n', ...
        'Scenario', 'd(eps_inf)', 'd(Delta_eps)', 'd(f_relax)', 'd(d_mm)');
    fprintf('%s\n', repmat('-', 1, 75));

    recovery = cell(3, 1);
    for iS = 1:3
        p0 = scenarios(iS, :);
        rec = zeros(nTrials, 4);
        for iT = 1:nTrials
            rng((iS-1)*100 + iT);
            T_gen  = debye_T(p0(1), p0(2), p0(3)*1e9, p0(4), f, c0);
            mag_n  = 20*log10(abs(T_gen)) + noise_mag_dB * randn(size(f));
            ph_n   = angle(T_gen) + deg2rad(noise_phase_deg) * randn(size(f));
            T_noisy = 10.^(mag_n/20) .* exp(1j * ph_n);

            [prec, ~] = debye_invert(T_noisy, f, lb_debye, ub_debye, ...
                                     popSize, maxIter, c0);
            rec(iT, :) = prec;
        end
        recovery{iS} = rec;
        err = mean(rec - p0, 1);
        fprintf('%-25s %+10.3f %+10.3f %+10.1f %+10.1f\n', ...
            scNames{iS}, err(1), err(2), (err(3))/1e9, err(4)*1000);
    end

    % ---- Plot recovery quality ----
    fig3 = figure('Color', 'w', 'Position', [50, 100, 1200, 400]);

    for iS = 1:3
        subplot(1, 3, iS);
        p0 = scenarios(iS, :);
        T_gen = debye_T(p0(1), p0(2), p0(3)*1e9, p0(4), f, c0);

        rng(42);
        mag_n = 20*log10(abs(T_gen)) + noise_mag_dB * randn(size(f));
        ph_n  = angle(T_gen) + deg2rad(noise_phase_deg) * randn(size(f));
        T_noisy = 10.^(mag_n/20) .* exp(1j * ph_n);

        [prec, ~] = debye_invert(T_noisy, f, lb_debye, ub_debye, ...
                                 popSize, maxIter, c0);
        T_fit = debye_T(prec(1), prec(2), prec(3), prec(4), f, c0);

        plot(fGHz, mag_n, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 4); hold on;
        plot(fGHz, 20*log10(abs(T_gen)), 'k-', 'LineWidth', 1.8);
        plot(fGHz, 20*log10(abs(T_fit)), 'r--', 'LineWidth', 1.5);

        xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 12);
        ylabel('|S_{21}| (dB)', 'FontName', fontName, 'FontSize', 12);
        title(sprintf('%s\nTrue: ei=%.1f De=%.1f fr=%.0fG d=%.0fmm\nFit:  ei=%.1f De=%.1f fr=%.0fG d=%.0fmm', ...
            scNames{iS}, p0(1), p0(2), p0(3), p0(4)*1000, ...
            prec(1), prec(2), prec(3)/1e9, prec(4)*1000), ...
            'FontName', fontName, 'FontSize', 9);
        set(gca, 'FontName', fontName, 'FontSize', 10, 'Box', 'on');
        legend({'Noisy', 'Truth', 'Debye fit'}, ...
            'FontName', fontName, 'FontSize', 8, 'Location', 'southwest');
        grid on;
    end

    sgtitle('Debye Model: Recovery Test with Low Noise (0.1 dB, 1 deg)', ...
        'FontName', fontName, 'FontSize', 14, 'FontWeight', 'bold');
    exportgraphics(fig3, fullfile(outDir, 'Exp4_Debye_Recovery.png'), 'Resolution', 300);
    fprintf('  Saved: Exp4_Debye_Recovery.png\n');

    % =================================================================
    % PART 4: Debye vs Power-Law Head-to-Head
    % =================================================================
    fprintf('Part 4: Debye vs Power-law comparison...\n');

    % Generate Debye data, fit with BOTH models
    pD = [2.5, 2.0, 15e9, 0.120];  % Debye truth
    T_debye_true = debye_T(pD(1), pD(2), pD(3), pD(4), f, c0);
    mag_dB = 20*log10(abs(T_debye_true));
    ph = angle(T_debye_true);
    mag_n = mag_dB + 0.1 * randn(size(f));
    T_noisy = 10.^(mag_n/20) .* exp(1j * ph);

    % Fit with Debye model
    [pD_fit, lossD] = debye_invert(T_noisy, f, lb_debye, ub_debye, ...
                                   40, 300, c0);
    T_D_fit = debye_T(pD_fit(1), pD_fit(2), pD_fit(3), pD_fit(4), f, c0);

    % Fit with power-law model
    lb_pl = [1.5, 1e-3, -1.5, 0.05];
    ub_pl = [8.0, 0.10,  1.5, 0.40];
    [pPL_fit, lossPL] = pl_invert(T_noisy, f, lb_pl, ub_pl, 40, 300, c0, 30e9);
    T_PL_fit = fabry_perot_T_pl(pPL_fit(1), pPL_fit(2), pPL_fit(3), ...
                                pPL_fit(4), f, c0, 30e9);

    fig4 = figure('Color', 'w', 'Position', [50, 100, 1100, 700]);

    % S21 magnitude fit
    subplot(2, 2, 1);
    plot(fGHz, 20*log10(abs(T_debye_true)), 'k-', 'LineWidth', 1.8); hold on;
    plot(fGHz, 20*log10(abs(T_D_fit)), 'b--', 'LineWidth', 1.5);
    plot(fGHz, 20*log10(abs(T_PL_fit)), 'r-.', 'LineWidth', 1.5);
    xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 12);
    ylabel('|S_{21}| (dB)', 'FontName', fontName, 'FontSize', 12);
    legend({'Debye truth', 'Debye fit', 'Power-law fit'}, ...
        'FontName', fontName, 'FontSize', 10, 'Location', 'southwest');
    title('|S_{21}| Magnitude Fit Comparison', ...
        'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', 11, 'Box', 'on'); grid on;

    % Residual
    subplot(2, 2, 2);
    res_D  = 20*log10(abs(T_D_fit))  - 20*log10(abs(T_debye_true));
    res_PL = 20*log10(abs(T_PL_fit)) - 20*log10(abs(T_debye_true));
    plot(fGHz, res_D, 'b-', 'LineWidth', 1.5); hold on;
    plot(fGHz, res_PL, 'r-', 'LineWidth', 1.5);
    yline(0, 'k--');
    xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 12);
    ylabel('Residual (dB)', 'FontName', fontName, 'FontSize', 12);
    legend({sprintf('Debye (RMS=%.3f dB)', rms(res_D)), ...
            sprintf('Power-law (RMS=%.3f dB)', rms(res_PL))}, ...
        'FontName', fontName, 'FontSize', 10, 'Location', 'best');
    title('Fitting Residual Comparison', ...
        'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', 11, 'Box', 'on'); grid on;

    % Convergence
    subplot(2, 2, 3);
    semilogy(lossD, 'b-', 'LineWidth', 1.5); hold on;
    semilogy(lossPL, 'r-', 'LineWidth', 1.5);
    xlabel('Iteration', 'FontName', fontName, 'FontSize', 12);
    ylabel('Cost', 'FontName', fontName, 'FontSize', 12);
    legend({sprintf('Debye (final=%.4f)', lossD(end)), ...
            sprintf('Power-law (final=%.4f)', lossPL(end))}, ...
        'FontName', fontName, 'FontSize', 10, 'Location', 'northeast');
    title('Convergence Comparison', ...
        'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', 11, 'Box', 'on'); grid on;

    % Effective permittivity comparison
    subplot(2, 2, 4);
    epsD_debye  = debye_eps(pD(1), pD(2), pD(3), f);
    epsD_fit_debye = debye_eps(pD_fit(1), pD_fit(2), pD_fit(3), f);
    epsPL_debye = fabry_perot_eps_pl(pPL_fit(1), pPL_fit(2), pPL_fit(3), f, 30e9);

    yyaxis left;
    plot(fGHz, real(epsD_debye), 'k-', 'LineWidth', 1.8); hold on;
    plot(fGHz, real(epsD_fit_debye), 'b--', 'LineWidth', 1.5);
    plot(fGHz, real(epsPL_debye), 'r-.', 'LineWidth', 1.5);
    ylabel('\epsilon'' (real part)', 'FontName', fontName, 'FontSize', 12);

    yyaxis right;
    plot(fGHz, -imag(epsD_debye), 'k-', 'LineWidth', 1.8); hold on;
    plot(fGHz, -imag(epsD_fit_debye), 'b--', 'LineWidth', 1.5);
    plot(fGHz, -imag(epsPL_debye), 'r-.', 'LineWidth', 1.5);
    ylabel('\epsilon'''' (imag part)', 'FontName', fontName, 'FontSize', 12);

    xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 12);
    title('Recovered \epsilon^*(f) Comparison', ...
        'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', 11, 'Box', 'on'); grid on;

    sgtitle('Debye vs Power-Law: Head-to-Head Model Comparison', ...
        'FontName', fontName, 'FontSize', 15, 'FontWeight', 'bold');
    exportgraphics(fig4, fullfile(outDir, 'Exp4_Debye_vs_PowerLaw.png'), 'Resolution', 300);
    fprintf('  Saved: Exp4_Debye_vs_PowerLaw.png\n');

    % ---- Summary table ----
    fprintf('\n=== Comparison Summary ===\n');
    fprintf('%-18s %8s %8s %10s %8s\n', 'Parameter', 'Truth', 'Debye', 'PowerLaw', '');
    fprintf('%-18s %8s %8s %10s %8s\n', '----------', '-----', '------', '---------', '');
    fprintf('%-18s %8.2f %8.2f %10.2f %8s\n', 'eps_inf / eps_r', pD(1), pD_fit(1), pPL_fit(1), '');
    fprintf('%-18s %8.2f %8.2f %10s %8s\n', 'Delta_eps', pD(2), pD_fit(2), '-- (N/A)', '');
    fprintf('%-18s %8.0f %8.0f %10s %8s\n', 'f_relax / n', pD(3)/1e9, pD_fit(3)/1e9, ...
        sprintf('n=%.2f', pPL_fit(3)), '');
    fprintf('%-18s %8.0f %8.0f %10.0f %8s\n', 'd (mm)', pD(4)*1000, pD_fit(4)*1000, pPL_fit(4)*1000, '');
    fprintf('%-18s %8s %8.3f %10.3f %8s\n', '|S21| RMS (dB)', '', ...
        rms(20*log10(abs(T_D_fit)) - mag_dB), ...
        rms(20*log10(abs(T_PL_fit)) - mag_dB), '');

    % ---- Save all data ----
    save(fullfile(outDir, 'Exp4_Debye_Data.mat'), ...
        'scenarios', 'scNames', 'recovery', 'pD', 'pD_fit', 'pPL_fit', ...
        'lossD', 'lossPL', 'curv');
    fprintf('\n=== Exp4 complete. Results in: %s ===\n', outDir);
end

%% ==================== DEBYE MODEL FUNCTIONS ====================

function eps_c = debye_eps(eps_inf, delta_eps, f_relax, f)
    eps_c = eps_inf + delta_eps ./ (1 + 1j * f / f_relax);
end

function T = debye_T(eps_inf, delta_eps, f_relax, d, f, c0)
    eps_c = debye_eps(eps_inf, delta_eps, f_relax, f);
    gamma = 1j * (2*pi*f) / c0 .* sqrt(eps_c);
    R = (sqrt(eps_c) - 1) ./ (sqrt(eps_c) + 1);
    E2 = exp(-2 * gamma * d);
    T = exp(-gamma * d) .* (1 - R.^2) ./ (1 - R.^2 .* E2);
end

%% ==================== POWER-LAW MODEL FUNCTIONS ====================

function T = fabry_perot_T_pl(eps_r, tand_ref, n, d, f, c0, fref)
    tand_f = tand_ref .* (f ./ fref).^n;
    eps_c  = eps_r .* (1 - 1j .* tand_f);
    gamma  = 1j * (2*pi*f) / c0 .* sqrt(eps_c);
    R = (sqrt(eps_c) - 1) ./ (sqrt(eps_c) + 1);
    E2 = exp(-2 * gamma * d);
    T = exp(-gamma * d) .* (1 - R.^2) ./ (1 - R.^2 .* E2);
end

function eps_c = fabry_perot_eps_pl(eps_r, tand_ref, n, f, fref)
    tand_f = tand_ref .* (f ./ fref).^n;
    eps_c  = eps_r .* (1 - 1j .* tand_f);
end

%% ==================== 1D SWEEP ====================

function cost = sweep_debye_1d(name, vals, p, f, c0, mag_true_dB)
    N = length(vals);
    cost = zeros(N, 1);
    for i = 1:N
        pp = p;
        switch name
            case 'eps_inf',   pp.eps_inf   = vals(i);
            case 'delta_eps', pp.delta_eps = vals(i);
            case 'f_relax',   pp.f_relax   = vals(i);
            case 'd',         pp.d         = vals(i);
        end
        T = debye_T(pp.eps_inf, pp.delta_eps, pp.f_relax, pp.d, f, c0);
        cost(i) = sqrt(mean((20*log10(abs(T)) - mag_true_dB).^2));
    end
end

%% ==================== 2D SWEEP ====================

function [C, mv1, mv2] = sweep_debye_2d(name1, vals1, name2, vals2, p, f, c0, mag_true_dB)
    N1 = length(vals1); N2 = length(vals2);
    C = zeros(N1, N2);
    for i = 1:N1
        for j = 1:N2
            pp = p;
            switch name1
                case 'eps_inf',   pp.eps_inf   = vals1(i);
                case 'delta_eps', pp.delta_eps = vals1(i);
                case 'f_relax',   pp.f_relax   = vals1(i);
                case 'd',         pp.d         = vals1(i);
            end
            switch name2
                case 'eps_inf',   pp.eps_inf   = vals2(j);
                case 'delta_eps', pp.delta_eps = vals2(j);
                case 'f_relax',   pp.f_relax   = vals2(j);
                case 'd',         pp.d         = vals2(j);
            end
            T = debye_T(pp.eps_inf, pp.delta_eps, pp.f_relax, pp.d, f, c0);
            C(i,j) = sqrt(mean((20*log10(abs(T)) - mag_true_dB).^2));
        end
    end
    [~, idx] = min(C(:));
    [ri, rj] = ind2sub(size(C), idx);
    mv1 = vals1(ri); mv2 = vals2(rj);
end

%% ==================== COUPLING RATIO ====================

function r = coupling_ratio(C, x, y)
    [minVal, idx] = min(C(:));
    [ri, rj] = ind2sub(size(C), idx);
    mask = C <= 2 * minVal;
    if ~any(mask(:)), r = nan; return; end
    [rows, cols] = find(mask);
    span_x = (max(x(rows)) - min(x(rows))) / (max(x) - min(x));
    span_y = (max(y(cols)) - min(y(cols))) / (max(y) - min(y));
    r = min(span_x, span_y) / max(span_x, span_y);
end

%% ==================== DEBYE INVERSION (DE) ====================

function [best_x, loss_curve] = debye_invert(T_meas, f, lb, ub, popSize, maxIter, c0)
    D = 4;
    mag_meas_dB = 20 * log10(abs(T_meas));
    phase_meas  = angle(T_meas);

    X = repmat(lb, popSize, 1) + lhsdesign(popSize, D) .* repmat(ub - lb, popSize, 1);

    cost = zeros(popSize, 1);
    for i = 1:popSize
        cost(i) = cost_debye(X(i,:), f, mag_meas_dB, phase_meas, c0);
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
            cost_u = cost_debye(u, f, mag_meas_dB, phase_meas, c0);
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

function J = cost_debye(x, f, mag_meas_dB, phase_meas, c0)
    T = debye_T(x(1), x(2), x(3), x(4), f, c0);
    mag_err   = 20*log10(abs(T)) - mag_meas_dB;
    phase_err = angle(T) - phase_meas;
    phase_err = wrapToPi(phase_err);
    L_mag   = mag_err.^2;   maskM = abs(mag_err) > 1.0;
    L_mag(maskM) = 2.0 * abs(mag_err(maskM)) - 1.0;
    L_phase = phase_err.^2; maskP = abs(phase_err) > 0.5;
    L_phase(maskP) = 2.0 * abs(phase_err(maskP)) - 0.25;
    J = 0.5 * mean(L_mag) + 0.25 * mean(L_phase);
end

%% ==================== POWER-LAW INVERSION (DE) ====================

function [best_x, loss_curve] = pl_invert(T_meas, f, lb, ub, popSize, maxIter, c0, fref)
    D = 4;
    mag_meas_dB = 20 * log10(abs(T_meas));
    phase_meas  = angle(T_meas);
    X = repmat(lb, popSize, 1) + lhsdesign(popSize, D) .* repmat(ub - lb, popSize, 1);

    cost = zeros(popSize, 1);
    for i = 1:popSize
        cost(i) = cost_pl(X(i,:), f, mag_meas_dB, phase_meas, c0, fref);
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
            cost_u = cost_pl(u, f, mag_meas_dB, phase_meas, c0, fref);
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

function J = cost_pl(x, f, mag_meas_dB, phase_meas, c0, fref)
    T = fabry_perot_T_pl(x(1), x(2), x(3), x(4), f, c0, fref);
    mag_err   = 20*log10(abs(T)) - mag_meas_dB;
    phase_err = angle(T) - phase_meas;
    phase_err = wrapToPi(phase_err);
    L_mag   = mag_err.^2;   maskM = abs(mag_err) > 1.0;
    L_mag(maskM) = 2.0 * abs(mag_err(maskM)) - 1.0;
    L_phase = phase_err.^2; maskP = abs(phase_err) > 0.5;
    L_phase(maskP) = 2.0 * abs(phase_err(maskP)) - 0.25;
    J = 0.5 * mean(L_mag) + 0.25 * mean(L_phase);
end
