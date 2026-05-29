%% Layered Wood Dielectric Inversion
% Transfer Matrix Method (TMM) with Chebyshev acceleration.
% Based on the full derivation in tmm_derivation.md.
%
% Model: alternating earlywood (E) / latewood (L) layers.
% Debye relaxation: eps*(f) = eps_inf + delta_eps / (1 + j*f/f_relax)
%
% 3-parameter optimization (with physically-motivated reductions):
%   X = [eps_inf_E, eps_inf_L, f_relax]
% Fixed: d_E, d_L (annual ring measurement), delta_eps_E, delta_eps_L
%
% Two polarization modes measured separately:
%   Parallel:  E-field along fiber (L-axis) -> probes eps_L*
%   Perpendicular: E-field across fiber (T-axis) -> probes eps_T*

function run_tmm_inversion()
    clc; close all;

    outDir = fullfile(pwd, 'simulations', 'results');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    fontName = 'Times New Roman';

    % ====== Physical Constants ======
    c0 = 299792458;
    Z0 = 377;

    % ====== Frequency Setup ======
    fGHz = linspace(18, 44, 201)';
    f    = fGHz * 1e9;

    % ====== Fixed Geometric Parameters ======
    % From caliper measurement: total thickness ~120mm (small wood)
    % From annual ring observation: ~3mm earlywood + ~2mm latewood per pair
    d_total  = 0.120;       % 120 mm total
    d_E_mm   = 3.0;         % earlywood layer thickness (mm)
    d_L_mm   = 2.0;         % latewood layer thickness (mm)
    d_E      = d_E_mm * 1e-3;
    d_L      = d_L_mm * 1e-3;
    N_pairs  = round(d_total / (d_E + d_L));
    N_layers = 2 * N_pairs;

    fprintf('=== Layered Wood TMM Inversion ===\n');
    fprintf('Total thickness: %.0f mm\n', d_total * 1000);
    fprintf('Layer pair: %.1f mm (E) + %.1f mm (L) = %.1f mm\n', ...
        d_E_mm, d_L_mm, d_E_mm + d_L_mm);
    fprintf('N_pairs = %d,  N_layers = %d\n\n', N_pairs, N_layers);

    % ====== Fixed Debye Parameters (from literature) ======
    % Earlywood: higher moisture -> larger delta_eps
    % Latewood:  lower moisture  -> smaller delta_eps
    %
    % Shared f_relax justification:
    %   Microscopically, earlywood and latewood differ in pore structure
    %   and bound/free water ratio. However, in the 18-44 GHz band,
    %   the dominant polarization mechanism is the strong Debye relaxation
    %   of free water. Bound water contribution in this band is largely
    %   saturated and absorbed into the high-frequency background eps_inf.
    %   Therefore, we assign a single equivalent f_relax representing the
    %   average relaxation frequency of free water in the wood matrix.
    %   This approximation preserves the physics while reducing the
    %   optimization dimension by 1, improving convergence and uniqueness.
    delta_eps_E = 3.0;
    delta_eps_L = 1.5;

    % ====== Load Measurement Data ======
    dataDir = fullfile(pwd, 'lzmwoods');

    % Air reference
    airFiles = dir(fullfile(dataDir, 'air*.csv'));
    airRef = load_measurement(fullfile(airFiles(1).folder, airFiles(1).name));
    S21_air = airRef.S_complex;

    % Wood sample files
    matFiles = dir(fullfile(dataDir, '*.csv'));
    matFiles = matFiles(~contains({matFiles.name}, 'air', 'IgnoreCase', true));

    % ====== Optimizer Settings ======
    % 3 parameters: [eps_inf_E, eps_inf_L, f_relax]
    lb3 = [1.5, 1.5, 0.5e9];
    ub3 = [10.0, 10.0, 60e9];
    popSize = 60;
    maxIter = 400;

    % ====== Process Each Sample ======
    fprintf('%-25s %10s %10s %12s %8s\n', ...
        'Sample', 'eps_inf_E', 'eps_inf_L', 'f_relax(GHz)', 'RMS(dB)');
    fprintf('%s\n', repmat('-', 1, 75));

    results = struct();

    for iF = 1:length(matFiles)
        fname = matFiles(iF).name;
        [~, baseName] = fileparts(fname);

        % ---- Load and preprocess ----
        mat = load_measurement(fullfile(matFiles(iF).folder, fname));

        % Interpolate to common frequency grid
        f_interp = fGHz * 1e9;
        S_mat = interp1(mat.fHz, mat.S_complex, f_interp, 'linear', 0);
        S_air = interp1(airRef.fHz, S21_air, f_interp, 'linear', 0);
        S_norm = S_mat ./ S_air;

        % Truncate: keep band where |S21| > -80 dB
        mask = 20*log10(abs(S_norm)) > -80;
        if sum(mask) < 50
            warning('Too few valid points for %s, skipping.', fname);
            continue;
        end
        f_use = f_interp(mask);
        S_use = S_norm(mask);

        % Downsample for optimization speed (~300 points)
        step = max(1, floor(length(f_use) / 300));
        f_fit = f_use(1:step:end);
        S_fit = S_use(1:step:end);

        % Savitzky-Golay smoothing (magnitude only)
        wlen = max(3, 2*floor(length(S_fit)/16) + 1);
        S_mag_sm = sgolayfilt(abs(S_fit), 2, wlen);
        S_fit_sm = S_mag_sm .* exp(1j * angle(S_fit));

        % ---- DE Optimization (3-param TMM + Debye) ----
        rng(iF * 100);
        [p_best, loss_curve] = tmm_debye_invert(...
            S_fit_sm, f_fit, lb3, ub3, ...
            delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, ...
            popSize, maxIter, c0, Z0);

        eps_inf_E = p_best(1);
        eps_inf_L = p_best(2);
        f_relax   = p_best(3);

        % Compute best-fit S21
        S21_model = tmm_debye_S21(eps_inf_E, eps_inf_L, f_relax, ...
            delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, f_fit, c0, Z0);
        rms_val = sqrt(mean((20*log10(abs(S21_model)) - 20*log10(abs(S_fit_sm))).^2));

        % ---- Store ----
        results(iF).name       = baseName;
        results(iF).eps_inf_E  = eps_inf_E;
        results(iF).eps_inf_L  = eps_inf_L;
        results(iF).f_relax    = f_relax;
        results(iF).rms        = rms_val;
        results(iF).f_fit      = f_fit;
        results(iF).S_meas     = S_fit_sm;
        results(iF).S_model    = S21_model;
        results(iF).loss_curve = loss_curve;

        fprintf('%-25s %10.2f %10.2f %12.1f %8.3f\n', ...
            baseName, eps_inf_E, eps_inf_L, f_relax/1e9, rms_val);
    end

    % =================================================================
    % FIGURE 1: S21 Fits
    % =================================================================
    fig1 = figure('Color', 'w', 'Position', [50, 100, 1100, 750]);

    for iF = 1:length(results)
        subplot(2, 2, iF);
        rr = results(iF);
        fplot = rr.f_fit / 1e9;

        plot(fplot, 20*log10(abs(rr.S_meas)), '.', ...
            'Color', [0.3 0.3 0.3], 'MarkerSize', 4); hold on;
        plot(fplot, 20*log10(abs(rr.S_model)), '-', ...
            'Color', [0.85 0.17 0.15], 'LineWidth', 1.8);

        xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 12);
        ylabel('|S_{21}| (dB)', 'FontName', fontName, 'FontSize', 12);
        title(sprintf(['%s (RMS=%.3f dB)\n\\epsilon_\\infty^E=%.2f, ' ...
            '\\epsilon_\\infty^L=%.2f, f_{relax}=%.1f GHz'], ...
            rr.name, rr.rms, rr.eps_inf_E, rr.eps_inf_L, rr.f_relax/1e9), ...
            'FontName', fontName, 'FontSize', 10);
        set(gca, 'FontName', fontName, 'FontSize', 11, 'Box', 'on');
        legend({'Measured', 'TMM+Debye'}, ...
            'FontName', fontName, 'FontSize', 9, 'Location', 'southwest');
        grid on;
    end

    sgtitle('Layered Wood TMM + Debye Inversion (3-parameter)', ...
        'FontName', fontName, 'FontSize', 14, 'FontWeight', 'bold');
    exportgraphics(fig1, fullfile(outDir, 'TMM_Inversion_Fits.png'), 'Resolution', 300);
    fprintf('\n  Saved: TMM_Inversion_Fits.png\n');

    % =================================================================
    % FIGURE 2: Texture Comparison
    % =================================================================
    fig2 = figure('Color', 'w', 'Position', [50, 100, 1000, 420]);

    pairs_cfg = {{'0_垂直', '0_平行'}, {'1_垂直', '1_平行'}};
    pair_titles = {'Large Wood', 'Small Wood'};
    colors  = {[0.29 0.49 0.73], [0.86 0.65 0.23]};

    for iPair = 1:2
        subplot(1, 2, iPair); hold on;

        for iS = 1:2
            target = pairs_cfg{iPair}{iS};
            for iF = 1:length(results)
                if strcmp(results(iF).name, target)
                    rr = results(iF);
                    % Plot measured
                    plot(rr.f_fit/1e9, 20*log10(abs(rr.S_meas)), '.', ...
                        'Color', colors{iS}, 'MarkerSize', 2);
                    % Plot model
                    plot(rr.f_fit/1e9, 20*log10(abs(rr.S_model)), '-', ...
                        'Color', colors{iS}, 'LineWidth', 1.8);
                    break;
                end
            end
        end

        xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 12);
        ylabel('|S_{21}| (dB)', 'FontName', fontName, 'FontSize', 12);
        title(pair_titles{iPair}, 'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');
        legend({'Vertical (meas)', 'Vertical (model)', ...
                'Parallel (meas)', 'Parallel (model)'}, ...
            'FontName', fontName, 'FontSize', 8, 'Location', 'southwest');
        set(gca, 'FontName', fontName, 'FontSize', 11, 'Box', 'on'); grid on;
    end

    sgtitle('Texture Effect: TMM Layered Model Fits', ...
        'FontName', fontName, 'FontSize', 14, 'FontWeight', 'bold');
    exportgraphics(fig2, fullfile(outDir, 'TMM_Texture_Comparison.png'), 'Resolution', 300);
    fprintf('  Saved: TMM_Texture_Comparison.png\n');

    % =================================================================
    % Summary Table
    % =================================================================
    fprintf('\n=============== Inversion Summary ===============\n');
    fprintf('%-22s %10s %10s %10s %8s\n', ...
        'Sample', 'eps_inf_E', 'eps_inf_L', 'f_relax(GHz)', 'RMS(dB)');
    fprintf('%s\n', repmat('-', 1, 70));
    for iF = 1:length(results)
        rr = results(iF);
        fprintf('%-22s %10.2f %10.2f %10.1f %8.3f\n', ...
            rr.name, rr.eps_inf_E, rr.eps_inf_L, rr.f_relax/1e9, rr.rms);
    end

    % Texture effect analysis
    fprintf('\n--- Texture Effect ---\n');
    for iPair = 1:2
        vert_name = pairs_cfg{iPair}{1};
        para_name = pairs_cfg{iPair}{2};
        vert_val = []; para_val = [];
        for iF = 1:length(results)
            if strcmp(results(iF).name, vert_name)
                vert_val = [results(iF).eps_inf_E, results(iF).eps_inf_L, results(iF).f_relax];
            end
            if strcmp(results(iF).name, para_name)
                para_val = [results(iF).eps_inf_E, results(iF).eps_inf_L, results(iF).f_relax];
            end
        end
        if ~isempty(vert_val) && ~isempty(para_val)
            fprintf('%s:\n', pair_titles{iPair});
            fprintf('  Vertical:  eps_inf_E=%.2f  eps_inf_L=%.2f  f_relax=%.1f GHz\n', ...
                vert_val(1), vert_val(2), vert_val(3)/1e9);
            fprintf('  Parallel:  eps_inf_E=%.2f  eps_inf_L=%.2f  f_relax=%.1f GHz\n', ...
                para_val(1), para_val(2), para_val(3)/1e9);
            fprintf('  Delta:     eps_inf_E=%.2f  eps_inf_L=%.2f  f_relax=%.1f GHz\n', ...
                vert_val(1)-para_val(1), vert_val(2)-para_val(2), ...
                (vert_val(3)-para_val(3))/1e9);
        end
    end

    save(fullfile(outDir, 'TMM_Inversion_Results.mat'), 'results', ...
        'd_total', 'd_E_mm', 'd_L_mm', 'N_pairs', 'delta_eps_E', 'delta_eps_L');
    fprintf('\nResults saved to: %s\n', fullfile(outDir, 'TMM_Inversion_Results.mat'));
end

%% ==================== TMM + DEBYE FORWARD MODEL ====================

function S21 = tmm_debye_S21(eps_inf_E, eps_inf_L, f_relax, ...
        delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, f, c0, Z0)
    % Compute S21 for alternating E-L layered structure.
    % Each layer uses Debye dispersion: eps*(f) = eps_inf + delta_eps/(1+j*f/f_relax)
    % Uses Chebyshev acceleration: O(1) per frequency, independent of N_pairs.

    eps_E = eps_inf_E + delta_eps_E ./ (1 + 1j * f / f_relax);
    eps_L = eps_inf_L + delta_eps_L ./ (1 + 1j * f / f_relax);

    [A, B, C, D] = tmm_chebyshev(eps_E, eps_L, d_E, d_L, N_pairs, f, c0, Z0);

    denom = A + B/Z0 + C*Z0 + D;
    S21 = 2 ./ denom;
end

function [A, B, C, D] = tmm_chebyshev(eps_E, eps_L, d_E, d_L, N_pairs, f, c0, Z0)
    % Compute total ABCD matrix using Chebyshev polynomial identity.
    % All inputs are N_freq × 1 vectors (eps_E, eps_L may be complex).
    % d_E, d_L are scalars.
    % Returns A, B, C, D as N_freq × 1 column vectors.

    n_E = sqrt(eps_E);  % complex refractive index, earlywood
    n_L = sqrt(eps_L);  % complex refractive index, latewood

    k0 = 2 * pi * f / c0;

    % Electrical lengths
    delta_E = k0 .* n_E * d_E;
    delta_L = k0 .* n_L * d_L;

    % Wave impedances
    Z_E = Z0 ./ n_E;
    Z_L = Z0 ./ n_L;

    % Layer characteristic matrices
    cE = cos(delta_E); sE = sin(delta_E);
    cL = cos(delta_L); sL = sin(delta_L);

    % M_pair = M_E * M_L, computed element-wise for all frequencies
    alpha = cE.*cL - (Z_E./Z_L).*sE.*sL;
    beta  = 1j * (Z_E.*sE.*cL + Z_L.*cE.*sL);
    gamma = 1j * ((1./Z_E).*sE.*cL + (1./Z_L).*cE.*sL);
    delta_term = cE.*cL - (Z_L./Z_E).*sE.*sL;
    % Note: det(M_pair) = alpha.*delta_term - beta.*gamma = 1 (verified)

    % Half-trace: s = (alpha + delta)/2
    s = (alpha + delta_term) / 2;

    % Chebyshev: U_{N-1}(s) and U_{N-2}(s) via theta = acos(s)
    theta = acos(s);

    % Guard against 0/0 = NaN when |s| -> 1 (sin(theta) -> 0)
    % Limit: U_n(1) = n+1, U_n(-1) = (-1)^n * (n+1)
    % For lossy wood, s is complex; this guards the rare real-s boundary case
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

    % M_total = U_{N-1} * M_pair - U_{N-2} * I
    A = U_Nm1 .* alpha      - U_Nm2;
    B = U_Nm1 .* beta;
    C = U_Nm1 .* gamma;
    D = U_Nm1 .* delta_term - U_Nm2;
end

%% ==================== DE OPTIMIZER (3-PARAM) ====================

function [best_x, loss_curve] = tmm_debye_invert(S_meas, f, lb, ub, ...
        delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, ...
        popSize, maxIter, c0, Z0)
    % Differential Evolution for 3-parameter TMM+Debye model.
    % X = [eps_inf_E, eps_inf_L, f_relax]

    D = 3;
    mag_meas_dB = 20 * log10(max(abs(S_meas), 1e-12));

    % Latin Hypercube initialization
    X = repmat(lb, popSize, 1) + lhsdesign(popSize, D) .* repmat(ub - lb, popSize, 1);

    cost = zeros(popSize, 1);
    for i = 1:popSize
        cost(i) = cost_tmm3(X(i,:), f, mag_meas_dB, ...
            delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, c0, Z0);
    end

    [best_cost, idx] = min(cost);
    best_x = X(idx, :);
    loss_curve = zeros(maxIter, 1);

    F  = 0.5;
    CR = 0.7;

    for iter = 1:maxIter
        for i = 1:popSize
            % DE/rand/1 mutation
            candidates = setdiff(1:popSize, i);
            r = candidates(randperm(length(candidates), 3));
            v = X(r(1), :) + F * (X(r(2), :) - X(r(3), :));

            % Binomial crossover
            j_rand = randi(D);
            u = X(i, :);
            for j = 1:D
                if rand < CR || j == j_rand
                    u(j) = v(j);
                end
            end

            % Boundary reflection
            u = max(lb, min(ub, u));

            % Selection
            cost_u = cost_tmm3(u, f, mag_meas_dB, ...
                delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, c0, Z0);
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

function J = cost_tmm3(x, f, mag_meas_dB, ...
        delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, c0, Z0)
    % Cost function: Huber loss on |S21| magnitude

    S21 = tmm_debye_S21(x(1), x(2), x(3), ...
        delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, f, c0, Z0);

    mag_err = 20*log10(max(abs(S21), 1e-12)) - mag_meas_dB;

    % Huber loss: quadratic for small errors, linear for large
    L = mag_err.^2;
    mask = abs(mag_err) > 1.0;
    L(mask) = 2.0 * abs(mag_err(mask)) - 1.0;

    J = mean(L);
end

%% ==================== CSV DATA LOADER ====================

function out = load_measurement(fp)
    assert(exist(fp,'file')==2, 'File not found: %s', fp);
    L = local_read_lines(fp);
    if isempty(L), error('Empty: %s', fp); end
    L{1} = local_strip_bom(L{1});

    keyRegex = '(Freq|Frequency).*(Hz)|S21\(DB\)|S21_DB|S21\(DEG\)|S21_DEG|Phase';
    idxH = local_find_contains(L, 'BEGIN');
    if ~isempty(idxH)
        c = idxH+1;
        if c <= numel(L) && ~isempty(regexp(L{c}, keyRegex, 'once', 'ignorecase'))
            idxH = c;
        else
            i2 = local_find_regex(L, keyRegex, idxH+1);
            if ~isempty(i2), idxH = i2; end
        end
    else
        idxH = local_find_regex(L, keyRegex, 1);
    end
    assert(~isempty(idxH), 'Header not found: %s', fp);

    hdr = strtrim(L{idxH});
    delim = local_guess_delim(hdr);
    names = local_split_header(hdr, delim);
    ncol = numel(names);

    dlines = L(idxH+1:end);
    mask = true(size(dlines));
    for i = 1:numel(dlines)
        s = strtrim(dlines{i});
        if isempty(s) || (~isempty(s) && s(1)=='!'), mask(i) = false; end
    end
    dlines = dlines(mask);

    cols = cell(1,ncol);
    for j = 1:ncol, cols{j} = strings(0,1); end
    for i = 1:numel(dlines)
        row = strsplit(dlines{i}, delim);
        if numel(row)>=ncol, row=row(1:ncol); else, row(end+1:ncol)={''}; end
        for j = 1:ncol
            cols{j}(end+1,1) = local_strip_quotes(strtrim(string(row{j})));
        end
    end

    wF = {'Freq_Hz_','Freq(Hz)','Freq_Hz','Freq','Frequency'};
    wA = {'S21(DB)','S21_DB','S21 dB','S21dB'};
    wP = {'S21(DEG)','S21_DEG','Phase','S21 DEG'};

    iF = pick_col(names, wF); iA = pick_col(names, wA); iP = pick_col(names, wP);

    fHz = str2double_clean(cols{iF});
    aDB = str2double_clean(cols{iA});
    pDeg = str2double_clean(cols{iP});

    mag = 10.^(aDB(:)/20);
    ph  = deg2rad(pDeg(:));
    S   = mag .* exp(1j*ph);
    [fHz, S] = local_clean(fHz(:), S(:));
    out = struct('fHz', fHz, 'S_complex', S);
end

function L = local_read_lines(fp)
    fid = fopen(fp,'r','n','UTF-8');
    if fid<0, error('Cannot open: %s', fp); end
    C = {}; i = 0;
    while true
        t = fgetl(fid); if ~ischar(t), break; end
        i=i+1; C{i,1}=t;
    end
    fclose(fid);
    L = C;
end

function s = local_strip_bom(s)
    if isempty(s), return; end
    u8 = uint8(s);
    if numel(u8)>=3 && isequal(u8(1:3),uint8([239 187 191])), s=char(u8(4:end)); end
    if ~isempty(s) && s(1)==char(65279), s=s(2:end); end
end

function idx = local_find_contains(L, token)
    idx = [];
    for i = 1:numel(L)
        if contains(L{i}, token, 'IgnoreCase', true), idx = i; return; end
    end
end

function idx = local_find_regex(L, re, st)
    if nargin<3, st=1; end
    idx = [];
    for i = st:numel(L)
        if ~isempty(regexp(L{i}, re, 'once', 'ignorecase')), idx = i; return; end
    end
end

function d = local_guess_delim(hdr)
    if contains(hdr, sprintf('\t')), d = sprintf('\t'); return; end
    c = [sum(hdr==','), sum(hdr==sprintf('\t')), sum(hdr==';')];
    [~,k] = max(c);
    if k==1, d=','; elseif k==2, d=sprintf('\t'); else, d=';'; end
end

function names = local_split_header(hdr, delim)
    P = strsplit(hdr, delim);
    names = cellfun(@(s) char(local_strip_quotes(string(strtrim(s)))), P, 'UniformOutput', false);
end

function s = local_strip_quotes(s)
    if strlength(s)>=2
        if extractBetween(s,1,1)=='"' && extractBetween(s,strlength(s),strlength(s))=='"'
            s = extractBetween(s,2,strlength(s)-1);
        end
    end
end

function v = str2double_clean(col)
    if iscell(col), col = string(col); end
    col = strrep(col,",",""); col = strrep(col,"，","");
    col = replace(col,'"','');
    v = str2double(col);
end

function idx = pick_col(names, cand)
    idx = [];
    for i = 1:numel(cand)
        j = find(strcmpi(names, cand{i}), 1, 'first');
        if ~isempty(j), idx = j; return; end
    end
    if isempty(idx), error('Column not found: %s', strjoin(cand,', ')); end
end

function [x, y] = local_clean(x, y)
    m = isfinite(x) & isfinite(real(y)) & isfinite(imag(y)) & x>0;
    x = x(m); y = y(m);
    if isempty(x), x=1; y=0; return; end
    [x, idx] = sort(x); y = y(idx);
    [ux,~,ic] = unique(x);
    if numel(ux) < numel(x)
        yr = accumarray(ic, [real(y), imag(y)], [], @(v) mean(v,1,'omitnan'));
        y = complex(yr(:,1), yr(:,2)); x = ux;
    end
    if numel(x)==1, x = x+[-1;1]*max(1,x*1e-9); y=y([1 1]); end
end
