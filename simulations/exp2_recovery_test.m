%% Experiment 2: Synthetic Data Recovery Test
% Generate S21 with known parameters, add noise, run inversion,
% compare recovered vs true values to assess algorithm accuracy.

function exp2_recovery_test()
    clc;

    outDir = fullfile(pwd, 'simulations', 'results');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    fontName = 'Times New Roman';

    c0   = 299792458;
    fref = 30e9;
    fGHz = linspace(18, 44, 201)';
    f    = fGHz * 1e9;

    % ---- 4 wood-type parameter sets ----
    % [eps_r, tand_ref, n, d(m)]
    truth = [
        3.50,  0.020, -0.50, 0.300;   % LargeWood_Vertical
        3.20,  0.025, -0.40, 0.300;   % LargeWood_Parallel
        4.00,  0.035, -0.60, 0.120;   % SmallWood_Vertical
        3.80,  0.040, -0.50, 0.120;   % SmallWood_Parallel
    ];
    sampleNames = {'LargeWood\_Vertical', 'LargeWood\_Parallel', ...
                   'SmallWood\_Vertical', 'SmallWood\_Parallel'};
    nSamples = size(truth, 1);

    % ---- Noise levels ----
    noiseLevels = {...
        struct('name','Noise-Free',   'mag_dB', 0.0,  'phase_deg', 0.0); ...
        struct('name','Low Noise',    'mag_dB', 0.1,  'phase_deg', 1.0); ...
        struct('name','High Noise',   'mag_dB', 0.5,  'phase_deg', 3.0); ...
    };
    nNoise = length(noiseLevels);
    nTrials = 5;  % Monte Carlo trials per condition

    % ---- Inversion bounds ----
    lb = [1.5,  1e-3, -1.5, 0.05];
    ub = [8.0,  0.10,  1.5, 0.50];

    % ---- DE inversion settings (compact) ----
    popSize  = 30;
    maxIter  = 200;
    F  = 0.5;
    CR = 0.7;

    % ---- Storage ----
    recovery_all = zeros(nSamples, nNoise, nTrials, 4);   % [samp, noise, trial, param]

    fprintf('=== Exp2: Synthetic Recovery Test ===\n');
    fprintf('%-22s %-14s %5s  %10s %10s %10s %10s\n', ...
        'Sample', 'Noise', 'Trial', 'eps_r_err', 'tand_err', 'n_err', 'd_err(mm)');
    fprintf('%s\n', repmat('-', 1, 90));

    for iS = 1:nSamples
        p_true = truth(iS, :);

        for iN = 1:nNoise
            nl = noiseLevels{iN};

            for iT = 1:nTrials
                % ---- Generate noisy synthetic data ----
                T_clean = fabry_perot_T(p_true(1), p_true(2), p_true(3), ...
                                        p_true(4), f, c0, fref);
                mag_dB_clean = 20 * log10(abs(T_clean));
                phase_clean  = angle(T_clean);

                rng((iS-1)*100 + (iN-1)*10 + iT);  % reproducible
                mag_dB_noisy   = mag_dB_clean + nl.mag_dB * randn(size(mag_dB_clean));
                phase_noisy    = phase_clean  + deg2rad(nl.phase_deg) * randn(size(phase_clean));
                T_noisy = 10.^(mag_dB_noisy/20) .* exp(1j * phase_noisy);

                % ---- DE inversion ----
                [p_rec, ~] = de_invert(T_noisy, f, lb, ub, popSize, maxIter, F, CR, c0, fref);

                recovery_all(iS, iN, iT, :) = p_rec;

                if iT <= 3  % print first 3 trials
                    fprintf('%-22s %-14s %5d  %+10.3f %+10.4f %+10.2f %+10.1f\n', ...
                        sampleNames{iS}, nl.name, iT, ...
                        p_rec(1)-p_true(1), p_rec(2)-p_true(2), ...
                        p_rec(3)-p_true(3), (p_rec(4)-p_true(4))*1000);
                end
            end
        end
    end

    % ---- Summary statistics ----
    fprintf('\n=== Recovery Error Summary (mean +/- std over %d trials) ===\n', nTrials);
    fprintf('%-22s %-14s  %10s %10s %10s %10s\n', ...
        'Sample', 'Noise', 'Err(eps_r)', 'Err(tand)', 'Err(n)', 'Err(d mm)');
    fprintf('%s\n', repmat('-', 1, 90));

    for iS = 1:nSamples
        p0 = truth(iS, :);
        for iN = 1:nNoise
            nl = noiseLevels{iN};
            err = reshape(recovery_all(iS, iN, :, :), nTrials, 4) - p0;
            err(:,4) = err(:,4) * 1000;  % d in mm
            mu = mean(err, 1);
            sd = std(err, 0, 1);
            fprintf('%-22s %-14s  %6.3f+/-%.3f %7.4f+/-%.4f %6.2f+/-%.2f %6.1f+/-%.1f\n', ...
                sampleNames{iS}, nl.name, mu(1), sd(1), mu(2), sd(2), mu(3), sd(3), mu(4), sd(4));
        end
        fprintf('%s\n', repmat('-', 1, 90));
    end

    % ---- Plot: recovery error bar chart ----
    fig1 = figure('Color', 'w', 'Position', [50, 100, 1400, 320]);

    paramNames = {'\epsilon_r''', 'tan\delta_{ref}', 'n', 'd (mm)'};
    paramUnits = {'', '', '', 'mm'};
    colors = {[0.2 0.4 0.8], [0.2 0.7 0.3], [0.9 0.4 0.2]};

    for iP = 1:4
        subplot(1, 4, iP);
        hold on;

        xPos = 1;
        tickLabels = {};
        tickPos = [];

        for iN = 1:nNoise
            nl = noiseLevels{iN};
            err_all = [];
            for iS = 1:nSamples
                p0 = truth(iS, :);
                err = reshape(recovery_all(iS, iN, :, iP), nTrials, 1) - p0(iP);
                if iP == 4, err = err * 1000; end
                err_all = [err_all; err];
            end
            mu  = mean(err_all);
            sd  = std(err_all);

            bar(xPos, mu, 'FaceColor', colors{iN}, 'EdgeColor', 'k', 'LineWidth', 0.8);
            errorbar(xPos, mu, sd, 'k.', 'LineWidth', 1.2, 'CapSize', 8);

            tickLabels{end+1} = nl.name;
            tickPos(end+1) = xPos;
            xPos = xPos + 1;
        end

        set(gca, 'XTick', tickPos, 'XTickLabel', tickLabels, ...
            'FontName', fontName, 'FontSize', 10, 'Box', 'on');
        if iP == 2
            ylabel('Relative Error (%)', 'FontName', fontName, 'FontSize', 12);
        else
            ylabel('Absolute Error', 'FontName', fontName, 'FontSize', 12);
        end
        title(paramNames{iP}, 'FontName', fontName, 'FontSize', 13, 'FontWeight', 'bold');
        grid on;
        yline(0, 'k-', 'LineWidth', 0.5);
    end

    sgtitle('Parameter Recovery Error by Noise Level (pooled across 4 samples, 10 trials each)', ...
        'FontName', fontName, 'FontSize', 14, 'FontWeight', 'bold');

    exportgraphics(fig1, fullfile(outDir, 'Exp2_Recovery_Error.png'), 'Resolution', 300);
    fprintf('\n  Saved: Exp2_Recovery_Error.png\n');

    % ---- Plot: S21 fitting examples ----
    fig2 = figure('Color', 'w', 'Position', [50, 100, 1300, 500]);

    % Pick one sample, show clean + noisy + fit
    iShow = 1;  % LargeWood_Vertical
    p_true = truth(iShow, :);
    T_clean = fabry_perot_T(p_true(1), p_true(2), p_true(3), p_true(4), f, c0, fref);
    mag_clean_dB = 20 * log10(abs(T_clean));

    for iN = 1:nNoise
        subplot(1, 3, iN);
        nl = noiseLevels{iN};

        rng(42);
        mag_noisy = mag_clean_dB + nl.mag_dB * randn(size(mag_clean_dB));
        phase_noisy = angle(T_clean) + deg2rad(nl.phase_deg) * randn(size(f));
        T_noisy = 10.^(mag_noisy/20) .* exp(1j * phase_noisy);

        [p_rec, ~] = de_invert(T_noisy, f, lb, ub, popSize, maxIter, F, CR, c0, fref);
        T_fit = fabry_perot_T(p_rec(1), p_rec(2), p_rec(3), p_rec(4), f, c0, fref);

        plot(fGHz, mag_noisy, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 5); hold on;
        plot(fGHz, mag_clean_dB, 'k-', 'LineWidth', 1.5);
        plot(fGHz, 20*log10(abs(T_fit)), 'r--', 'LineWidth', 1.5);

        xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 12);
        ylabel('|S_{21}| (dB)', 'FontName', fontName, 'FontSize', 12);
        title(sprintf('%s\nRecovered: er=%.2f tand=%.3f n=%.2f d=%.0fmm', ...
            nl.name, p_rec(1), p_rec(2), p_rec(3), p_rec(4)*1000), ...
            'FontName', fontName, 'FontSize', 11);
        set(gca, 'FontName', fontName, 'FontSize', 10, 'Box', 'on');
        legend({'Noisy data', 'True model', 'Recovered fit'}, ...
            'FontName', fontName, 'FontSize', 9, 'Location', 'southwest');
        grid on;
    end

    sgtitle(sprintf('Synthetic Recovery Example: %s', strrep(sampleNames{iShow}, '\', '')), ...
        'FontName', fontName, 'FontSize', 14, 'FontWeight', 'bold');

    exportgraphics(fig2, fullfile(outDir, 'Exp2_Fitting_Example.png'), 'Resolution', 300);
    fprintf('  Saved: Exp2_Fitting_Example.png\n');

    % ---- Save data ----
    save(fullfile(outDir, 'Exp2_Recovery_Data.mat'), ...
        'recovery_all', 'truth', 'sampleNames', 'noiseLevels', 'nTrials');
    fprintf('\n=== Exp2 complete. Results in: %s ===\n', outDir);
end

%% ====== Forward Model ======
function T = fabry_perot_T(eps_r, tand_ref, n, d, f, c0, fref)
    tand_f = tand_ref .* (f ./ fref).^n;
    eps_c  = eps_r .* (1 - 1j .* tand_f);
    omega  = 2 * pi * f;
    gamma  = 1j .* omega ./ c0 .* sqrt(eps_c);
    E = exp(-gamma .* d);
    R = (sqrt(eps_c) - 1) ./ (sqrt(eps_c) + 1);
    T = E .* (1 - R.^2) ./ (1 - R.^2 .* E.^2);
end

%% ====== Differential Evolution Inversion ======
function [best_x, loss_curve] = de_invert(T_meas, f, lb, ub, popSize, maxIter, F, CR, c0, fref)
    D = 4;
    mag_meas_dB = 20 * log10(abs(T_meas));
    phase_meas  = angle(T_meas);

    % LHS initialization
    X = repmat(lb, popSize, 1) + lhsdesign(popSize, D) .* repmat(ub - lb, popSize, 1);

    cost = zeros(popSize, 1);
    for i = 1:popSize
        cost(i) = compute_cost(X(i,:), f, mag_meas_dB, phase_meas, c0, fref);
    end

    [best_cost, idx] = min(cost);
    best_x = X(idx, :);
    loss_curve = zeros(maxIter, 1);

    for iter = 1:maxIter
        for i = 1:popSize
            % Mutation
            candidates = setdiff(1:popSize, i);
            r = candidates(randperm(length(candidates), 3));
            v = X(r(1), :) + F * (X(r(2), :) - X(r(3), :));

            % Crossover
            j_rand = randi(D);
            u = X(i, :);
            for j = 1:D
                if rand < CR || j == j_rand
                    u(j) = v(j);
                end
            end

            % Bound handling (reflection)
            u = max(lb, min(ub, u));

            % Selection
            cost_u = compute_cost(u, f, mag_meas_dB, phase_meas, c0, fref);
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

function J = compute_cost(x, f, mag_meas_dB, phase_meas, c0, fref)
    T = fabry_perot_T(x(1), x(2), x(3), x(4), f, c0, fref);
    mag_err  = 20*log10(abs(T)) - mag_meas_dB;
    phase_err = angle(T) - phase_meas;
    phase_err = wrapToPi(phase_err);

    % Huber-like loss
    L_mag = mag_err.^2;
    maskM = abs(mag_err) > 1.0;
    L_mag(maskM) = 2.0 * abs(mag_err(maskM)) - 1.0;
    L_phase = phase_err.^2;
    maskP = abs(phase_err) > 0.5;
    L_phase(maskP) = 2.0 * abs(phase_err(maskP)) - 0.25;

    J = 0.5 * mean(L_mag) + 0.25 * mean(L_phase);
end
