%% Experiment 5: Fixed-Thickness Debye Inversion
% With d independently measured, reduce from 4 to 3 parameters.
% Test on synthetic data first, then apply to real measurements.

function exp5_fixed_d()
    clc; close all;

    outDir = fullfile(pwd, 'simulations', 'results');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    fontName = 'Times New Roman';
    c0 = 299792458;
    fGHz = linspace(18, 44, 201)';
    f    = fGHz * 1e9;

    % =================================================================
    % PART 1: Synthetic Recovery with Fixed d (3-param Debye)
    % =================================================================
    fprintf('=== Part 1: Fixed-d Debye Recovery on Synthetic Data ===\n\n');

    % True parameters for moist wood (small block)
    p_true = struct('eps_inf', 2.5, 'delta_eps', 2.0, 'f_relax', 15e9, 'd', 0.120);

    % Fix d to true value, invert 3 params only
    lb3 = [1.5, 0.5,  1e9];
    ub3 = [5.0, 6.0, 80e9];
    popSize = 40;
    maxIter = 300;
    nTrials = 10;

    noise_mag_dB = 0.1;
    noise_phase_deg = 1.0;

    % Test with d_fixed at true value AND with d_fixed slightly wrong
    d_fixed_cases = [
        p_true.d;            % exact
        p_true.d * 1.02;     % +2% error (~2.4mm for 120mm)
        p_true.d * 1.05;     % +5% error (~6mm)
    ];
    case_names = {'d fixed = true (120.0 mm)', ...
                  'd fixed = true + 2% (122.4 mm)', ...
                  'd fixed = true + 5% (126.0 mm)'};

    fprintf('%-28s %10s %10s %10s\n', 'Case', 'd(eps_inf)', 'd(Delta_eps)', 'd(f_relax)');
    fprintf('%s\n', repmat('-', 1, 68));

    for iC = 1:length(d_fixed_cases)
        d_fixed = d_fixed_cases(iC);
        rec = zeros(nTrials, 3);
        for iT = 1:nTrials
            rng((iC-1)*100 + iT);
            T_gen  = debye_T(p_true.eps_inf, p_true.delta_eps, ...
                             p_true.f_relax, p_true.d, f, c0);
            mag_n  = 20*log10(abs(T_gen)) + noise_mag_dB * randn(size(f));
            ph_n   = angle(T_gen) + deg2rad(noise_phase_deg) * randn(size(f));
            T_noisy = 10.^(mag_n/20) .* exp(1j * ph_n);

            [prec, ~] = debye3_invert(T_noisy, f, lb3, ub3, d_fixed, ...
                                      popSize, maxIter, c0);
            rec(iT, :) = prec;
        end
        err = mean(rec, 1) - [p_true.eps_inf, p_true.delta_eps, p_true.f_relax];
        err_std = std(rec, 0, 1);
        fprintf('%-28s %+10.3f %+10.3f %+10.1f\n', ...
            case_names{iC}, err(1), err(2), err(3)/1e9);
        fprintf('  (std)                  %10.3f %10.3f %10.1f\n', ...
            err_std(1), err_std(2), err_std(3)/1e9);
    end

    % =================================================================
    % PART 2: Head-to-Head: Fixed-d Debye vs Fixed-d Power-Law
    % =================================================================
    fprintf('\n=== Part 2: Fixed-d Debye vs Fixed-d Power-law ===\n\n');

    % Generate Debye data
    T_true = debye_T(p_true.eps_inf, p_true.delta_eps, p_true.f_relax, ...
                     p_true.d, f, c0);
    rng(42);
    mag_noisy = 20*log10(abs(T_true)) + 0.1 * randn(size(f));
    ph_noisy  = angle(T_true) + deg2rad(1.0) * randn(size(f));
    T_noisy = 10.^(mag_noisy/20) .* exp(1j * ph_noisy);

    % Fixed-d Debye fit (3 params)
    [pD3, lossD3] = debye3_invert(T_noisy, f, lb3, ub3, p_true.d, ...
                                  50, 400, c0);
    T_D3 = debye_T(pD3(1), pD3(2), pD3(3), p_true.d, f, c0);
    rms_D3 = rms(20*log10(abs(T_D3)) - 20*log10(abs(T_true)));

    % Fixed-d Power-law fit (3 params)
    lb3_pl = [1.5, 1e-3, -1.5];
    ub3_pl = [8.0, 0.10,  1.5];
    [pPL3, lossPL3] = pl3_invert(T_noisy, f, lb3_pl, ub3_pl, p_true.d, ...
                                 50, 400, c0, 30e9);
    T_PL3 = fabry_perot_T_pl(pPL3(1), pPL3(2), pPL3(3), p_true.d, f, c0, 30e9);
    rms_PL3 = rms(20*log10(abs(T_PL3)) - 20*log10(abs(T_true)));

    fprintf('%-22s %8s %8s %10s %8s %8s\n', ...
        'Model', 'eps_inf/r', 'Delta/n', 'f_relax(GHz)', 'd(mm)', 'RMS(dB)');
    fprintf('%s\n', repmat('-', 1, 80));
    fprintf('%-22s %8s %8s %10s %8s %8s\n', ...
        'Truth (Debye)', '2.50', '2.00', '15.0', '120', '--');
    fprintf('%-22s %8.2f %8.2f %10.1f %8.0f %8.3f\n', ...
        'Debye fit (3-param)', pD3(1), pD3(2), pD3(3)/1e9, p_true.d*1000, rms_D3);
    fprintf('%-22s %8.2f %8s %10s %8.0f %8.3f\n', ...
        'PowerLaw fit (3-param)', pPL3(1), sprintf('n=%.2f',pPL3(3)), ...
        sprintf('tand=%.3f',pPL3(2)), p_true.d*1000, rms_PL3);

    % =================================================================
    % PART 3: Apply Fixed-d Debye to Real Measurement Data
    % =================================================================
    fprintf('\n=== Part 3: Fixed-d Debye on Real Wood Measurements ===\n\n');

    dataDir = fullfile(pwd, 'lzmwoods');
    airFiles = dir(fullfile(dataDir, 'air*.csv'));
    airRef = load_measurement(fullfile(airFiles(1).folder, airFiles(1).name));
    S21_air = airRef.S_complex;

    matFiles = dir(fullfile(dataDir, '*.csv'));
    matFiles = matFiles(~contains({matFiles.name}, 'air', 'IgnoreCase', true));

    % Sample thickness (user-measured)
    d_small  = 0.120;   % Small wood block: ~120mm
    d_large  = 0.300;   % Large wood block: ~300mm (approximate)

    fprintf('%-25s %8s %8s %10s %8s\n', ...
        'Sample', 'eps_inf', 'Delta_eps', 'f_relax(GHz)', 'RMS(dB)');
    fprintf('%s\n', repmat('-', 1, 75));

    results_real = struct();

    for iF = 1:length(matFiles)
        fname = matFiles(iF).name;
        [~, base, ~] = fileparts(fname);
        mat = load_measurement(fullfile(matFiles(iF).folder, fname));

        % Interpolate to air frequency grid
        f_interp = fGHz * 1e9;
        S_mat = interp1(mat.fHz, mat.S_complex, f_interp, 'linear', 0);
        S_air = interp1(airRef.fHz, S21_air, f_interp, 'linear', 0);
        S_norm = S_mat ./ S_air;

        % Truncate to valid band
        mask = 20*log10(abs(S_norm)) > -80;
        f_use  = f_interp(mask);
        S_use  = S_norm(mask);

        % Downsample for speed
        step = max(1, floor(length(f_use) / 300));
        f_fit = f_use(1:step:end);
        S_fit = S_use(1:step:end);

        % Smooth magnitude
        wlen = max(3, 2*floor(length(S_fit)/16) + 1);
        S_mag_sm = sgolayfilt(abs(S_fit), 2, wlen);
        S_fit_sm = S_mag_sm .* exp(1j * angle(S_fit));

        % Determine thickness
        if contains(base, '0_'), d_fixed = d_large;
        else, d_fixed = d_small; end

        % Fixed-d Debye inversion
        [pD, lossD] = debye3_invert(S_fit_sm, f_fit, lb3, ub3, d_fixed, ...
                                    50, 400, c0);
        T_fit = debye_T(pD(1), pD(2), pD(3), d_fixed, f_fit, c0);
        rms_val = rms(20*log10(abs(T_fit)) - 20*log10(abs(S_fit_sm)));

        results_real(iF).name    = base;
        results_real(iF).p       = pD;
        results_real(iF).d_fixed = d_fixed;
        results_real(iF).rms     = rms_val;
        results_real(iF).f_fit   = f_fit;
        results_real(iF).S_meas  = S_fit_sm;
        results_real(iF).T_fit   = T_fit;
        results_real(iF).loss    = lossD;

        fprintf('%-25s %8.2f %8.2f %10.1f %8.3f\n', ...
            base, pD(1), pD(2), pD(3)/1e9, rms_val);
    end

    % =================================================================
    % PART 4: Plot Real Data Fits
    % =================================================================
    fig1 = figure('Color', 'w', 'Position', [50, 100, 1100, 750]);

    for iF = 1:length(matFiles)
        subplot(2, 2, iF);
        rr = results_real(iF);
        fplot = rr.f_fit / 1e9;
        h1 = plot(fplot, 20*log10(abs(rr.S_meas)), '.', ...
            'Color', [0.3 0.3 0.3], 'MarkerSize', 5); hold on;
        h2 = plot(fplot, 20*log10(abs(rr.T_fit)), '-', ...
            'Color', [0.85 0.17 0.15], 'LineWidth', 1.8);

        xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 12);
        ylabel('|S_{21}| (dB)', 'FontName', fontName, 'FontSize', 12);
        p = rr.p;
        title(sprintf('%s (d=%.0fmm fixed)\nei=%.2f De=%.2f fr=%.1fG RMS=%.3fdB', ...
            rr.name, rr.d_fixed*1000, p(1), p(2), p(3)/1e9, rr.rms), ...
            'FontName', fontName, 'FontSize', 10);
        set(gca, 'FontName', fontName, 'FontSize', 11, 'Box', 'on');
        legend([h1, h2], {'Measured', 'Debye fit (3-param)'}, ...
            'FontName', fontName, 'FontSize', 9, 'Location', 'southwest');
        grid on;
    end

    sgtitle('Real Wood: Fixed-Thickness Debye Inversion (d from caliper)', ...
        'FontName', fontName, 'FontSize', 15, 'FontWeight', 'bold');
    exportgraphics(fig1, fullfile(outDir, 'Exp5_RealWood_FixedD_Debye.png'), ...
        'Resolution', 300);
    fprintf('\n  Saved: Exp5_RealWood_FixedD_Debye.png\n');

    % =================================================================
    % PART 5: Texture Comparison (Vertical vs Parallel)
    % =================================================================
    fig2 = figure('Color', 'w', 'Position', [50, 100, 1000, 420]);

    % Pair samples: large wood (0_) and small wood (1_)
    pairs_config = {{'0_垂直', '0_平行'}, {'1_垂直', '1_平行'}};
    pair_titles = {'Large Wood (d=300mm)', 'Small Wood (d=120mm)'};

    for iPair = 1:2
        subplot(1, 2, iPair);
        hold on;
        colors = {[0.29 0.49 0.73], [0.86 0.65 0.23]};  % blue, gold
        styles = {'-', '--'};
        labels = {'Vertical', 'Parallel'};
        legend_handles = [];

        for iS = 1:2
            target = pairs_config{iPair}{iS};
            for iF = 1:length(results_real)
                if strcmp(results_real(iF).name, target)
                    rr = results_real(iF);
                    h = plot(rr.f_fit/1e9, 20*log10(abs(rr.T_fit)), ...
                        styles{iS}, 'Color', colors{iS}, 'LineWidth', 1.8);
                    legend_handles = [legend_handles, h];
                    % Also plot measured
                    plot(rr.f_fit/1e9, 20*log10(abs(rr.S_meas)), '.', ...
                        'Color', colors{iS}, 'MarkerSize', 2);
                    break;
                end
            end
        end

        xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 12);
        ylabel('|S_{21}| (dB)', 'FontName', fontName, 'FontSize', 12);
        title(pair_titles{iPair}, 'FontName', fontName, 'FontSize', 13, ...
            'FontWeight', 'bold');
        set(gca, 'FontName', fontName, 'FontSize', 11, 'Box', 'on');
        legend(legend_handles, labels, 'FontName', fontName, ...
            'FontSize', 10, 'Location', 'southwest');
        grid on;
    end

    sgtitle('Texture Effect: Debye Fits for Vertical vs Parallel Grain', ...
        'FontName', fontName, 'FontSize', 15, 'FontWeight', 'bold');
    exportgraphics(fig2, fullfile(outDir, 'Exp5_Texture_Comparison.png'), ...
        'Resolution', 300);
    fprintf('  Saved: Exp5_Texture_Comparison.png\n');

    % ---- Save ----
    save(fullfile(outDir, 'Exp5_Results.mat'), 'results_real', 'd_small', 'd_large');
    fprintf('\n=== Exp5 complete. Results in: %s ===\n', outDir);
end

%% ==================== 3-PARAM DEBYE FUNCTIONS ====================

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

function [best_x, loss_curve] = debye3_invert(T_meas, f, lb, ub, d_fixed, ...
        popSize, maxIter, c0)
    D = 3;
    mag_meas_dB = 20 * log10(abs(T_meas));
    phase_meas  = angle(T_meas);

    X = repmat(lb, popSize, 1) + lhsdesign(popSize, D) .* repmat(ub - lb, popSize, 1);

    cost = zeros(popSize, 1);
    for i = 1:popSize
        cost(i) = cost3(X(i,:), f, mag_meas_dB, phase_meas, d_fixed, c0);
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
            cost_u = cost3(u, f, mag_meas_dB, phase_meas, d_fixed, c0);
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

function J = cost3(x, f, mag_meas_dB, phase_meas, d_fixed, c0)
    T = debye_T(x(1), x(2), x(3), d_fixed, f, c0);
    mag_err = 20*log10(abs(T)) - mag_meas_dB;
    ph_err  = angle(T) - phase_meas;
    ph_err  = wrapToPi(ph_err);
    Lm = mag_err.^2; mask = abs(mag_err) > 1.0;
    Lm(mask) = 2.0 * abs(mag_err(mask)) - 1.0;
    Lp = ph_err.^2; mask2 = abs(ph_err) > 0.5;
    Lp(mask2) = 2.0 * abs(ph_err(mask2)) - 0.25;
    J = 0.5 * mean(Lm) + 0.25 * mean(Lp);
end

%% ==================== 3-PARAM POWER-LAW FUNCTIONS ====================

function T = fabry_perot_T_pl(eps_r, tand_ref, n, d, f, c0, fref)
    tand_f = tand_ref .* (f ./ fref).^n;
    eps_c  = eps_r .* (1 - 1j .* tand_f);
    gamma  = 1j * (2*pi*f) / c0 .* sqrt(eps_c);
    R = (sqrt(eps_c) - 1) ./ (sqrt(eps_c) + 1);
    E2 = exp(-2 * gamma * d);
    T = exp(-gamma * d) .* (1 - R.^2) ./ (1 - R.^2 .* E2);
end

function [best_x, loss_curve] = pl3_invert(T_meas, f, lb, ub, d_fixed, ...
        popSize, maxIter, c0, fref)
    D = 3;
    mag_meas_dB = 20 * log10(abs(T_meas));
    phase_meas  = angle(T_meas);
    X = repmat(lb, popSize, 1) + lhsdesign(popSize, D) .* repmat(ub - lb, popSize, 1);
    cost = zeros(popSize, 1);
    for i = 1:popSize
        cost(i) = cost_pl3(X(i,:), f, mag_meas_dB, phase_meas, d_fixed, c0, fref);
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
            cost_u = cost_pl3(u, f, mag_meas_dB, phase_meas, d_fixed, c0, fref);
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

function J = cost_pl3(x, f, mag_meas_dB, phase_meas, d_fixed, c0, fref)
    T = fabry_perot_T_pl(x(1), x(2), x(3), d_fixed, f, c0, fref);
    mag_err = 20*log10(abs(T)) - mag_meas_dB;
    ph_err  = angle(T) - phase_meas;
    ph_err  = wrapToPi(ph_err);
    Lm = mag_err.^2; mask = abs(mag_err) > 1.0;
    Lm(mask) = 2.0 * abs(mag_err(mask)) - 1.0;
    Lp = ph_err.^2; mask2 = abs(ph_err) > 0.5;
    Lp(mask2) = 2.0 * abs(ph_err(mask2)) - 0.25;
    J = 0.5 * mean(Lm) + 0.25 * mean(Lp);
end

%% ==================== CSV LOADER (from project) ====================

function out = load_measurement(fp)
    assert(exist(fp,'file')==2, 'File not found: %s', fp);
    L = local_read_all_lines(fp);
    if isempty(L), error('Empty file: %s', fp); end
    L{1} = local_strip_bom(L{1});
    keyRegex = '(Freq|Frequency).*(Hz)|S21\(DB\)|S21_DB|S21\(DEG\)|S21_DEG|Phase';
    idxHeader = local_find_line_contains(L, 'BEGIN');
    if ~isempty(idxHeader)
        idxCand = idxHeader + 1;
        if idxCand <= numel(L) && ~isempty(regexp(L{idxCand}, keyRegex, 'once', 'ignorecase'))
            idxHeader = idxCand;
        else
            idx2 = local_find_line_regex(L, keyRegex, idxHeader+1);
            if ~isempty(idx2), idxHeader = idx2; end
        end
    else
        idxHeader = local_find_line_regex(L, keyRegex, 1);
    end
    assert(~isempty(idxHeader), 'Header not found: %s', fp);
    headerLine = strtrim(L{idxHeader});
    delim = local_guess_delim(headerLine);
    names = local_split_header(headerLine, delim);
    ncol  = numel(names);
    dataLines = L(idxHeader+1:end);
    mask = true(size(dataLines));
    for i=1:numel(dataLines)
        s = strtrim(dataLines{i});
        if isempty(s) || (~isempty(s) && s(1)=='!'), mask(i) = false; end
    end
    dataLines = dataLines(mask);
    cols = cell(1,ncol);
    for j=1:ncol, cols{j} = strings(0,1); end
    for i=1:numel(dataLines)
        row = strsplit(dataLines{i}, delim);
        if numel(row) >= ncol, row = row(1:ncol);
        else, row(end+1:ncol) = {''}; end
        for j=1:ncol
            cols{j}(end+1,1) = local_strip_quotes(strtrim(string(row{j})));
        end
    end
    wantF = {'Freq_Hz_','Freq(Hz)','Freq_Hz','Freq','Frequency','Frequency_Hz_','Frequency(Hz)'};
    wantA = {'S21(DB)','S21_DB_','S21_DB','S21 dB','S21dB'};
    wantP = {'S21(DEG)','S21_DEG_','S21_DEG','Phase','S21 DEG','S21deg'};
    iF = pick_alias_cs(names, wantF, false);
    if isempty(iF), iF = pick_alias_ci(names, wantF, 'Freq'); end
    iA = pick_alias_cs(names, wantA, false);
    if isempty(iA), iA = pick_alias_ci(names, wantA, 'S21 dB'); end
    iP = pick_alias_cs(names, wantP, false);
    if isempty(iP), iP = pick_alias_ci(names, wantP, 'S21 DEG'); end
    fHz  = str2double_clean(cols{iF});
    ADB  = str2double_clean(cols{iA});
    Pdeg = str2double_clean(cols{iP});
    mag = 10.^(ADB(:)/20); ph = deg2rad(Pdeg(:));
    S = mag .* exp(1j*ph);
    [fHz, S] = local_clean_series(fHz(:), S(:));
    out = struct('fHz', fHz, 'S_complex', S);
end

function L = local_read_all_lines(fp)
    fid = fopen(fp,'r','n','UTF-8');
    if fid<0, error('Cannot open: %s', fp); end
    C = {}; i = 0;
    while true
        t = fgetl(fid); if ~ischar(t), break; end
        i=i+1; C{i,1} = t;
    end
    fclose(fid);
    if isempty(C), L = {}; else, L = C; end
end
function s = local_strip_bom(s)
    if isempty(s), return; end
    u8 = uint8(s);
    if numel(u8)>=3 && isequal(u8(1:3), uint8([239 187 191])), s = char(u8(4:end)); end
    if ~isempty(s) && s(1)==char(65279), s = s(2:end); end
end
function idx = local_find_line_contains(L, token)
    idx = [];
    for i=1:numel(L)
        if contains(L{i}, token, 'IgnoreCase', true), idx = i; return; end
    end
end
function idx = local_find_line_regex(L, re, startIdx)
    if nargin<3, startIdx=1; end
    idx = [];
    for i=startIdx:numel(L)
        if ~isempty(regexp(L{i}, re, 'once', 'ignorecase')), idx = i; return; end
    end
end
function d = local_guess_delim(headerLine)
    if contains(headerLine, sprintf('\t')), d = sprintf('\t'); return; end
    c_comma = sum(headerLine==','); c_tab = sum(headerLine==sprintf('\t'));
    c_sc = sum(headerLine==';');
    [~,k] = max([c_comma,c_tab,c_sc]); d = ',';
    if k==2, d=sprintf('\t'); end; if k==3, d=';'; end
end
function names = local_split_header(headerLine, delim)
    P = strsplit(headerLine, delim);
    names = cellfun(@(s) char(local_strip_quotes(string(strtrim(s)))), P, 'UniformOutput', false);
end
function s = local_strip_quotes(s)
    if strlength(s)>=2
        if s.extractBetween(1,1) == '"' && s.extractBetween(strlength(s),strlength(s)) == '"'
            s = extractBetween(s, 2, strlength(s)-1);
        end
    end
end
function v = str2double_clean(col)
    if iscell(col), col = string(col); end
    col = strrep(col, ",", ""); col = strrep(col, "，", "");
    col = replace(col, '"', ""); v = str2double(col);
end
function idx = pick_alias_cs(names, cand, must)
    idx = [];
    for i=1:numel(cand)
        j = find(strcmp(names, cand{i}), 1, 'first');
        if ~isempty(j), idx = j; return; end
    end
    if nargin>=3 && must, error('Column not found: %s', strjoin(cand, ', ')); end
end
function idx = pick_alias_ci(names, cand, tag)
    idx = [];
    for i=1:numel(cand)
        j = find(strcmpi(names, cand{i}), 1, 'first');
        if ~isempty(j)
            fprintf('[Hint] %s matched case-insensitive: %s\n', tag, names{j});
            idx = j; return;
        end
    end
    error('Column not found (%s): %s', tag, strjoin(cand, ', '));
end
function [x, y] = local_clean_series(x, y)
    m = isfinite(x) & isfinite(real(y)) & isfinite(imag(y)) & x>0;
    x = x(m); y = y(m);
    if isempty(x), x = 1; y = 0; return; end
    [x, idx] = sort(x); y = y(idx);
    [ux, ~, ic] = unique(x);
    if numel(ux) < numel(x)
        yr = accumarray(ic, [real(y), imag(y)], [], @(v) mean(v,1,'omitnan'));
        y = complex(yr(:,1), yr(:,2)); x = ux;
    end
    if numel(x)==1
        x = x + [-1;1]*max(1,x*1e-9); y = y([1 1]);
    end
end
