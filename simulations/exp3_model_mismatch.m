%% Experiment 3: Model Mismatch Test
% Generate S21 using a more realistic Debye relaxation model,
% then invert with the simple power-law loss model.
% Quantifies the error introduced by using the wrong physics.

function exp3_model_mismatch()
    clc;

    outDir = fullfile(pwd, 'simulations', 'results');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    fontName = 'Times New Roman';
    c0   = 299792458;
    fref = 30e9;
    fGHz = linspace(18, 44, 201)';
    f    = fGHz * 1e9;

    % ---- Debye model parameters (more physical for moist wood) ----
    % eps*(f) = eps_inf + delta_eps / (1 + j * f / f_relax)
    % Then convert to effective (eps_r, tand) for comparison

    debye_true = [
        % [eps_inf, delta_eps, f_relax(GHz), d(m)]
        2.5,  2.0,  15.0, 0.300;   % Sample 1: relaxation at 15 GHz
        2.5,  1.5,  18.0, 0.300;   % Sample 2
        2.8,  2.5,  12.0, 0.120;   % Sample 3
        2.8,  2.0,  14.0, 0.120;   % Sample 4
    ];
    sampleNames = {'Sample1 (d=300mm)', 'Sample2 (d=300mm)', ...
                   'Sample3 (d=120mm)', 'Sample4 (d=120mm)'};
    nSamples = size(debye_true, 1);

    % ---- Noise level for this experiment (moderate) ----
    noise_mag_dB  = 0.1;
    noise_phase   = deg2rad(1.0);

    % ---- Power-law inversion bounds ----
    lb = [1.5,  1e-3, -1.5, 0.05];
    ub = [8.0,  0.10,  1.5, 0.50];
    popSize = 40;
    maxIter = 300;
    F  = 0.5;
    CR = 0.7;

    fprintf('=== Exp3: Model Mismatch Test ===\n');
    fprintf('Generating data with Debye relaxation, inverting with power-law model.\n\n');

    results = struct();
    nTrials = 10;

    for iS = 1:nSamples
        pD = debye_true(iS, :);
        eps_inf  = pD(1);
        delta_eps = pD(2);
        f_relax  = pD(3) * 1e9;
        d_true   = pD(4);

        % ---- Generate Debye S21 ----
        T_debye = debye_transmission(eps_inf, delta_eps, f_relax, d_true, f, c0);
        mag_debye_dB = 20 * log10(abs(T_debye));

        % ---- Effective permittivity and tand at reference ----
        eps_debye = eps_inf + delta_eps ./ (1 + 1j * f / f_relax);
        eps_at_ref = eps_inf + delta_eps / (1 + 1j * 30e9 / f_relax);
        eps_r_eff = real(eps_at_ref);
        tand_eff  = -imag(eps_at_ref) / real(eps_at_ref);

        fprintf('Sample %d: Debye -> eps_eff@30GHz = %.2f - j%.4f (tand=%.4f)\n', ...
            iS, eps_r_eff, eps_r_eff * tand_eff, tand_eff);

        % ---- Monte Carlo with noise ----
        p_rec_all = zeros(nTrials, 4);
        cost_all = zeros(nTrials, 1);

        for iT = 1:nTrials
            rng((iS-1)*100 + iT);
            mag_noisy  = mag_debye_dB + noise_mag_dB * randn(size(mag_debye_dB));
            phase_noisy = angle(T_debye) + noise_phase * randn(size(f));
            T_noisy = 10.^(mag_noisy/20) .* exp(1j * phase_noisy);

            [p_rec, loss] = de_invert(T_noisy, f, lb, ub, popSize, maxIter, F, CR, c0, fref);
            p_rec_all(iT, :) = p_rec;
            cost_all(iT) = loss(end);
        end

        results(iS).name    = sampleNames{iS};
        results(iS).debye   = pD;
        results(iS).eps_r_eff = eps_r_eff;
        results(iS).tand_eff  = tand_eff;
        results(iS).p_rec   = p_rec_all;
        results(iS).cost    = cost_all;
    end

    % ---- Plot: Power-law fit vs Debye truth ----
    fig1 = figure('Color', 'w', 'Position', [50, 100, 1300, 850]);

    for iS = 1:nSamples
        subplot(2, 2, iS);
        pD = debye_true(iS, :);
        T_debye = debye_transmission(pD(1), pD(2), pD(3)*1e9, pD(4), f, c0);

        % Best power-law fit (lowest cost)
        [~, iBest] = min(results(iS).cost);
        p_best = results(iS).p_rec(iBest, :);
        T_fit = fabry_perot_T(p_best(1), p_best(2), p_best(3), p_best(4), f, c0, fref);

        plot(fGHz, 20*log10(abs(T_debye)), 'k-', 'LineWidth', 2.0); hold on;
        plot(fGHz, 20*log10(abs(T_fit)), 'r--', 'LineWidth', 1.8);

        xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 12);
        ylabel('|S_{21}| (dB)', 'FontName', fontName, 'FontSize', 12);
        title(sprintf('%s', results(iS).name), ...
            'FontName', fontName, 'FontSize', 12, 'FontWeight', 'bold');
        set(gca, 'FontName', fontName, 'FontSize', 10, 'Box', 'on');
        legend({'Debye (truth)', 'Power-law fit'}, ...
            'FontName', fontName, 'FontSize', 9, 'Location', 'southwest');
        grid on;

        % Annotation
        txt = sprintf(['Debye: \\epsilon_\\infty=%.1f, \\Delta\\epsilon=%.1f, ' ...
            'f_{relax}=%.0fGHz\nPower-law fit: \\epsilon_r=%.2f, ' ...
            'tan\\delta=%.3f, n=%.2f\nd_{true}=%.0fmm, d_{fit}=%.0fmm'], ...
            pD(1), pD(2), pD(3), p_best(1), p_best(2), p_best(3), ...
            pD(4)*1000, p_best(4)*1000);
        text(0.98, 0.15, txt, 'Units', 'normalized', ...
            'FontName', fontName, 'FontSize', 8, ...
            'HorizontalAlignment', 'right', ...
            'BackgroundColor', [1 1 1 0.8], 'EdgeColor', [0.3 0.3 0.3]);
    end

    sgtitle('Model Mismatch: Debye-generated data fitted with Power-Law model', ...
        'FontName', fontName, 'FontSize', 14, 'FontWeight', 'bold');
    exportgraphics(fig1, fullfile(outDir, 'Exp3_Model_Mismatch_Fits.png'), 'Resolution', 300);
    fprintf('\n  Saved: Exp3_Model_Mismatch_Fits.png\n');

    % ---- Plot: Mismatch residual ----
    fig2 = figure('Color', 'w', 'Position', [50, 100, 1300, 400]);

    for iS = 1:nSamples
        subplot(1, 4, iS);
        pD = debye_true(iS, :);
        T_debye = debye_transmission(pD(1), pD(2), pD(3)*1e9, pD(4), f, c0);
        mag_debye = 20*log10(abs(T_debye));

        [~, iBest] = min(results(iS).cost);
        p_best = results(iS).p_rec(iBest, :);
        T_fit = fabry_perot_T(p_best(1), p_best(2), p_best(3), p_best(4), f, c0, fref);
        mag_fit = 20*log10(abs(T_fit));

        residual = mag_fit - mag_debye;
        plot(fGHz, residual, 'b-', 'LineWidth', 1.5); hold on;
        yline(0, 'k--');
        xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 11);
        ylabel('Residual (dB)', 'FontName', fontName, 'FontSize', 11);
        title(results(iS).name, 'FontName', fontName, 'FontSize', 11, 'FontWeight', 'bold');
        set(gca, 'FontName', fontName, 'FontSize', 10, 'Box', 'on'); grid on;
        rms_res = rms(residual);
        text(0.98, 0.92, sprintf('RMS=%.3f dB', rms_res), 'Units', 'normalized', ...
            'FontName', fontName, 'FontSize', 9, ...
            'HorizontalAlignment', 'right');
    end

    sgtitle('Mismatch Residual: Power-law fit minus Debye truth', ...
        'FontName', fontName, 'FontSize', 14, 'FontWeight', 'bold');
    exportgraphics(fig2, fullfile(outDir, 'Exp3_Mismatch_Residual.png'), 'Resolution', 300);
    fprintf('  Saved: Exp3_Mismatch_Residual.png\n');

    % ---- Print summary ----
    fprintf('\n=== Mismatch Summary ===\n');
    fprintf('%-25s %10s %10s %10s %10s\n', 'Sample', 'eps_r(fit)', 'tand(fit)', 'n(fit)', 'd(fit,mm)');
    for iS = 1:nSamples
        [~, iBest] = min(results(iS).cost);
        p = results(iS).p_rec(iBest, :);
        fprintf('%-25s %10.2f %10.3f %10.2f %10.0f\n', ...
            results(iS).name, p(1), p(2), p(3), p(4)*1000);
        fprintf('  True (effective):  %10.2f %10.3f  %10s %10.0f\n', ...
            results(iS).eps_r_eff, results(iS).tand_eff, '--', debye_true(iS,4)*1000);
    end

    save(fullfile(outDir, 'Exp3_Model_Mismatch_Data.mat'), 'results', 'debye_true');
    fprintf('\n=== Exp3 complete. Results in: %s ===\n', outDir);
end

%% ====== Debye Forward Model ======
function T = debye_transmission(eps_inf, delta_eps, f_relax, d, f, c0)
    eps_c = eps_inf + delta_eps ./ (1 + 1j * f / f_relax);
    omega = 2 * pi * f;
    gamma = 1j .* omega ./ c0 .* sqrt(eps_c);
    E = exp(-gamma .* d);
    R = (sqrt(eps_c) - 1) ./ (sqrt(eps_c) + 1);
    T = E .* (1 - R.^2) ./ (1 - R.^2 .* E.^2);
end

%% ====== Power-Law Forward Model ======
function T = fabry_perot_T(eps_r, tand_ref, n, d, f, c0, fref)
    tand_f = tand_ref .* (f ./ fref).^n;
    eps_c  = eps_r .* (1 - 1j .* tand_f);
    omega  = 2 * pi * f;
    gamma  = 1j .* omega ./ c0 .* sqrt(eps_c);
    E = exp(-gamma .* d);
    R = (sqrt(eps_c) - 1) ./ (sqrt(eps_c) + 1);
    T = E .* (1 - R.^2) ./ (1 - R.^2 .* E.^2);
end

%% ====== DE Inversion (reuse from exp2) ======
function [best_x, loss_curve] = de_invert(T_meas, f, lb, ub, popSize, maxIter, F, CR, c0, fref)
    D = 4;
    mag_meas_dB = 20 * log10(abs(T_meas));
    phase_meas  = angle(T_meas);
    X = repmat(lb, popSize, 1) + lhsdesign(popSize, D) .* repmat(ub - lb, popSize, 1);

    cost = zeros(popSize, 1);
    for i = 1:popSize
        cost(i) = mismatch_cost(X(i,:), f, mag_meas_dB, phase_meas, c0, fref);
    end

    [best_cost, idx] = min(cost);
    best_x = X(idx, :);
    loss_curve = zeros(maxIter, 1);

    for iter = 1:maxIter
        for i = 1:popSize
            candidates = setdiff(1:popSize, i);
            r = candidates(randperm(length(candidates), 3));
            v = X(r(1), :) + F * (X(r(2), :) - X(r(3), :));
            j_rand = randi(D);
            u = X(i, :);
            for j = 1:D
                if rand < CR || j == j_rand
                    u(j) = v(j);
                end
            end
            u = max(lb, min(ub, u));
            cost_u = mismatch_cost(u, f, mag_meas_dB, phase_meas, c0, fref);
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

function J = mismatch_cost(x, f, mag_meas_dB, phase_meas, c0, fref)
    T = fabry_perot_T(x(1), x(2), x(3), x(4), f, c0, fref);
    mag_err   = 20*log10(abs(T)) - mag_meas_dB;
    phase_err = angle(T) - phase_meas;
    phase_err = wrapToPi(phase_err);
    L_mag   = mag_err.^2;
    L_phase = phase_err.^2;
    maskM = abs(mag_err) > 1.0;
    L_mag(maskM) = 2.0 * abs(mag_err(maskM)) - 1.0;
    maskP = abs(phase_err) > 0.5;
    L_phase(maskP) = 2.0 * abs(phase_err(maskP)) - 0.25;
    J = 0.5 * mean(L_mag) + 0.25 * mean(L_phase);
end
