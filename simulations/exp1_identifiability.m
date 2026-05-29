%% Experiment 1: Parameter Identifiability Analysis
% Sweep each parameter to visualize cost landscape shape.
% Deep narrow V = identifiable; shallow basin = ill-conditioned.

function exp1_identifiability()
    clc; close all;

    outDir = fullfile(pwd, 'simulations', 'results');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    % ---- Reference parameters (representative wood values) ----
    p_true.eps_r     = 3.5;
    p_true.tand_ref  = 0.02;     % tan(delta) at 30 GHz
    p_true.n         = -0.5;
    p_true.d         = 0.120;    % 120 mm

    % ---- Frequency setup (matching measurement: 18-44 GHz) ----
    fGHz = linspace(18, 44, 201)';
    f    = fGHz * 1e9;
    c0   = 299792458;
    fref = 30e9;

    % ---- Generate synthetic "measured" S21 ----
    T_true = fabry_perot_T(p_true.eps_r, p_true.tand_ref, p_true.n, ...
                           p_true.d, f, c0, fref);
    mag_true_dB = 20 * log10(abs(T_true));

    % ---- Sweep ranges (physical bounds) ----
    eps_vals     = linspace(1.5, 8.0, 80);
    logtand_vals = linspace(-3.0, -0.1, 80);  % log10 scale
    n_vals       = linspace(-1.5, 1.5, 80);
    d_vals       = linspace(0.05, 0.50, 80);   % 50-500 mm

    tand_vals = 10.^logtand_vals;

    % =================================================================
    % 1D PARAMETER SWEEPS
    % =================================================================
    fprintf('Running 1D parameter sweeps...\n');

    cost_eps  = sweep_1d('eps_r',  eps_vals,     p_true, f, c0, fref, mag_true_dB);
    cost_tand = sweep_1d('tand_ref', tand_vals,  p_true, f, c0, fref, mag_true_dB);
    cost_n    = sweep_1d('n',        n_vals,     p_true, f, c0, fref, mag_true_dB);
    cost_d    = sweep_1d('d',        d_vals,     p_true, f, c0, fref, mag_true_dB);

    % ---- Plot 1D results ----
    fontName = 'Times New Roman';
    fontSizeAxes  = 12;
    fontSizeLabel = 14;
    fontSizeTitle = 14;

    fig1 = figure('Color', 'w', 'Position', [50, 100, 1000, 700]);

    subplot(2,2,1);
    plot(eps_vals, cost_eps, 'b-', 'LineWidth', 1.5); hold on;
    xline(p_true.eps_r, 'r--', 'LineWidth', 1.2);
    xlabel('\epsilon_r''', 'FontName', fontName, 'FontSize', fontSizeLabel);
    ylabel('RMSE (dB)', 'FontName', fontName, 'FontSize', fontSizeLabel);
    title('Cost vs \epsilon_r''', 'FontName', fontName, 'FontSize', fontSizeTitle, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', fontSizeAxes, 'Box', 'on'); grid on;

    subplot(2,2,2);
    semilogx(tand_vals, cost_tand, 'b-', 'LineWidth', 1.5); hold on;
    xline(p_true.tand_ref, 'r--', 'LineWidth', 1.2);
    xlabel('tan\delta @30GHz', 'FontName', fontName, 'FontSize', fontSizeLabel);
    ylabel('RMSE (dB)', 'FontName', fontName, 'FontSize', fontSizeLabel);
    title('Cost vs tan\delta', 'FontName', fontName, 'FontSize', fontSizeTitle, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', fontSizeAxes, 'Box', 'on'); grid on;

    subplot(2,2,3);
    plot(n_vals, cost_n, 'b-', 'LineWidth', 1.5); hold on;
    xline(p_true.n, 'r--', 'LineWidth', 1.2);
    xlabel('n (frequency exponent)', 'FontName', fontName, 'FontSize', fontSizeLabel);
    ylabel('RMSE (dB)', 'FontName', fontName, 'FontSize', fontSizeLabel);
    title('Cost vs n', 'FontName', fontName, 'FontSize', fontSizeTitle, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', fontSizeAxes, 'Box', 'on'); grid on;

    subplot(2,2,4);
    plot(d_vals*1000, cost_d, 'b-', 'LineWidth', 1.5); hold on;
    xline(p_true.d*1000, 'r--', 'LineWidth', 1.2);
    xlabel('d (mm)', 'FontName', fontName, 'FontSize', fontSizeLabel);
    ylabel('RMSE (dB)', 'FontName', fontName, 'FontSize', fontSizeLabel);
    title('Cost vs d', 'FontName', fontName, 'FontSize', fontSizeTitle, 'FontWeight', 'bold');
    set(gca, 'FontName', fontName, 'FontSize', fontSizeAxes, 'Box', 'on'); grid on;

    sgtitle('1D Parameter Sensitivity: Fabry-Perot Model (no noise)', ...
        'FontName', fontName, 'FontSize', 16, 'FontWeight', 'bold');

    exportgraphics(fig1, fullfile(outDir, 'Exp1_1D_Sensitivity.png'), 'Resolution', 300);
    fprintf('  Saved: Exp1_1D_Sensitivity.png\n');

    % ---- Curvature diagnostics ----
    fprintf('\n=== Curvature at true value (larger = better identified) ===\n');
    fprintf('eps_r:     %.4f\n', local_curvature(eps_vals, cost_eps, p_true.eps_r));
    fprintf('tan(delta):%.4f\n', local_curvature(tand_vals, cost_tand, p_true.tand_ref));
    fprintf('n:         %.4f\n', local_curvature(n_vals, cost_n, p_true.n));
    fprintf('d:         %.4f\n', local_curvature(d_vals, cost_d, p_true.d));

    % =================================================================
    % 2D COST CONTOURS: (eps_r, d) and (tand_ref, n)
    % =================================================================
    fprintf('\nRunning 2D contour sweeps...\n');

    % ---- (eps_r, d) contour ----
    [C_ed, minEps, minD] = sweep_2d('eps_r', eps_vals, 'd', d_vals, ...
                                     p_true, f, c0, fref, mag_true_dB);

    fig2 = figure('Color', 'w', 'Position', [50, 100, 1000, 420]);

    subplot(1,2,1);
    contourf(eps_vals, d_vals*1000, C_ed', 30, 'LineStyle', 'none');
    hold on;
    plot(p_true.eps_r, p_true.d*1000, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
    plot(minEps, minD*1000, 'w+', 'MarkerSize', 12, 'LineWidth', 2);
    xlabel('\epsilon_r''', 'FontName', fontName, 'FontSize', fontSizeLabel);
    ylabel('d (mm)', 'FontName', fontName, 'FontSize', fontSizeLabel);
    title('\epsilon_r'' vs d Coupling', 'FontName', fontName, ...
        'FontSize', fontSizeTitle, 'FontWeight', 'bold');
    colorbar; colormap(jet);
    legend({'','True value', 'Cost minimum'}, 'FontName', fontName, ...
        'FontSize', 10, 'Location', 'northeast');
    set(gca, 'FontName', fontName, 'FontSize', fontSizeAxes, 'Box', 'on');

    % ---- (tand_ref, n) contour ----
    [C_tn, minTand, minN] = sweep_2d('tand_ref', tand_vals, 'n', n_vals, ...
                                      p_true, f, c0, fref, mag_true_dB);

    subplot(1,2,2);
    contourf(log10(tand_vals), n_vals, C_tn', 30, 'LineStyle', 'none');
    hold on;
    plot(log10(p_true.tand_ref), p_true.n, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
    plot(log10(minTand), minN, 'w+', 'MarkerSize', 12, 'LineWidth', 2);
    xlabel('log_{10}(tan\delta_{ref})', 'FontName', fontName, 'FontSize', fontSizeLabel);
    ylabel('n', 'FontName', fontName, 'FontSize', fontSizeLabel);
    title('tan\delta_{ref} vs n Coupling', 'FontName', fontName, ...
        'FontSize', fontSizeTitle, 'FontWeight', 'bold');
    colorbar; colormap(jet);
    set(gca, 'FontName', fontName, 'FontSize', fontSizeAxes, 'Box', 'on');

    sgtitle('2D Cost Contours: Parameter Coupling Diagnostics', ...
        'FontName', fontName, 'FontSize', 16, 'FontWeight', 'bold');

    exportgraphics(fig2, fullfile(outDir, 'Exp1_2D_Coupling.png'), 'Resolution', 300);
    fprintf('  Saved: Exp1_2D_Coupling.png\n');

    % ---- Condition number estimate ----
    fprintf('\n=== Coupling strength (lower = stronger coupling) ===\n');
    fprintf('(eps_r, d) aspect ratio:       %.2f\n', ...
        coupling_aspect(C_ed, eps_vals, d_vals));
    fprintf('(tand_ref, n) aspect ratio:    %.2f\n', ...
        coupling_aspect(C_tn, log10(tand_vals), n_vals));

    fprintf('\n=== Exp1 complete. Results in: %s ===\n', outDir);
end

%% ====== Forward Model ======
function T = fabry_perot_T(eps_r, tand_ref, n, d, f, c0, fref)
    tand_f = tand_ref .* (f ./ fref).^n;
    eps_c  = eps_r .* (1 - 1j .* tand_f);

    omega = 2 * pi * f;
    gamma = 1j .* omega ./ c0 .* sqrt(eps_c);

    E = exp(-gamma .* d);
    R = (sqrt(eps_c) - 1) ./ (sqrt(eps_c) + 1);
    T = E .* (1 - R.^2) ./ (1 - R.^2 .* E.^2);
end

%% ====== 1D Sweep ======
function cost = sweep_1d(name, vals, p, f, c0, fref, mag_true_dB)
    N = length(vals);
    cost = zeros(N, 1);
    for i = 1:N
        pp = p;
        switch name
            case 'eps_r',     pp.eps_r    = vals(i);
            case 'tand_ref',  pp.tand_ref = vals(i);
            case 'n',         pp.n        = vals(i);
            case 'd',         pp.d        = vals(i);
        end
        T = fabry_perot_T(pp.eps_r, pp.tand_ref, pp.n, pp.d, f, c0, fref);
        cost(i) = sqrt(mean((20*log10(abs(T)) - mag_true_dB).^2));
    end
end

%% ====== 2D Sweep ======
function [C, minVal1, minVal2] = sweep_2d(name1, vals1, name2, vals2, p, f, c0, fref, mag_true_dB)
    N1 = length(vals1);
    N2 = length(vals2);
    C = zeros(N1, N2);
    for i = 1:N1
        for j = 1:N2
            pp = p;
            switch name1
                case 'eps_r',     pp.eps_r    = vals1(i);
                case 'tand_ref',  pp.tand_ref = vals1(i);
                case 'n',         pp.n        = vals1(i);
                case 'd',         pp.d        = vals1(i);
            end
            switch name2
                case 'eps_r',     pp.eps_r    = vals2(j);
                case 'tand_ref',  pp.tand_ref = vals2(j);
                case 'n',         pp.n        = vals2(j);
                case 'd',         pp.d        = vals2(j);
            end
            T = fabry_perot_T(pp.eps_r, pp.tand_ref, pp.n, pp.d, f, c0, fref);
            C(i,j) = sqrt(mean((20*log10(abs(T)) - mag_true_dB).^2));
        end
    end
    [minC, idx] = min(C(:));
    [ri, rj] = ind2sub(size(C), idx);
    minVal1 = vals1(ri);
    minVal2 = vals2(rj);
    fprintf('  2D min at (%.3f, %.3f), cost=%.4f\n', minVal1, minVal2, minC);
end

%% ====== Curvature (finite difference) ======
function k = local_curvature(x, y, x0)
    [~, i0] = min(abs(x - x0));
    i0 = max(2, min(i0, length(x)-1));
    dx = x(2) - x(1);
    dy = (y(i0+1) - 2*y(i0) + y(i0-1)) / dx^2;
    k  = abs(dy) / (1 + ((y(i0+1)-y(i0-1))/(2*dx))^2)^1.5;
end

%% ====== Coupling aspect ratio ======
function r = coupling_aspect(C, x, y)
    % Estimate the aspect ratio of the cost valley near the minimum.
    % r < 0.3 => strong coupling (elongated valley).
    % r > 0.7 => weak coupling (round bowl).
    [minVal, idx] = min(C(:));
    [ri, rj] = ind2sub(size(C), idx);

    % Collect half-width at 2x min cost
    mask = C <= 2 * minVal;
    if ~any(mask(:))
        r = nan; return;
    end
    [rows, cols] = find(mask);
    span_x = (max(x(rows)) - min(x(rows))) / (max(x) - min(x));
    span_y = (max(y(cols)) - min(y(cols))) / (max(y) - min(y));
    r = min(span_x, span_y) / max(span_x, span_y);
end
