function plot_posterior(results, outDir, fontName, ...
        delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, c0, Z0)
    % Generate 5 posterior diagnostic figures from MCMC results.

    n_samples = length(results);
    n_params = 3;
    param_names = {'\epsilon_\infty^E', '\epsilon_\infty^L', 'f_{relax} (GHz)'};
    n_chains = length(results(1).mcmc_chains);

    % ================================================================
    % FIGURE 1: Trace Plots (n_samples x 3 grid, 4 chains overlaid)
    % ================================================================
    fig1 = figure('Color', 'w', 'Position', [50, 50, 1100, 900]);
    chain_colors = lines(n_chains);

    for iF = 1:n_samples
        rr = results(iF);
        chains = rr.mcmc_chains;

        for iP = 1:n_params
            subplot(n_samples, n_params, (iF-1)*n_params + iP);
            hold on;
            for c = 1:n_chains
                y = chains{c}(:, iP);
                if iP == 3, y = y / 1e9; end
                plot(y, 'Color', [chain_colors(c,:) 0.4], 'LineWidth', 0.3);
            end
            ylabel(param_names{iP}, 'FontName', fontName, 'FontSize', 9);
            if iP == 1
                title(rr.name, 'FontName', fontName, 'FontSize', 10, ...
                    'FontWeight', 'bold', 'Interpreter', 'none');
            end
            set(gca, 'FontName', fontName, 'FontSize', 8, 'Box', 'on');
        end
    end
    sgtitle('MCMC Trace Plots', 'FontName', fontName, 'FontSize', 14, 'FontWeight', 'bold');
    exportgraphics(fig1, fullfile(outDir, 'Posterior_Trace.png'), 'Resolution', 150);
    fprintf('  Saved: Posterior_Trace.png\n');

    % ================================================================
    % FIGURE 2: Posterior Marginal Densities with 95% HDI
    % ================================================================
    fig2 = figure('Color', 'w', 'Position', [50, 50, 1100, 900]);

    for iF = 1:n_samples
        rr = results(iF);
        chains = rr.mcmc_chains;
        posterior = [];
        for c = 1:n_chains
            posterior = [posterior; chains{c}];
        end

        for iP = 1:n_params
            subplot(n_samples, n_params, (iF-1)*n_params + iP);
            hold on;

            vals = posterior(:, iP);
            if iP == 3, vals = vals / 1e9; end

            [f_est, xi] = ksdensity(vals);
            fill([xi, fliplr(xi)], [f_est, zeros(size(f_est))], ...
                [0.7 0.85 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.6);

            hdi_lo = quantile(vals, 0.025);
            hdi_hi = quantile(vals, 0.975);
            yl = ylim;
            plot([hdi_lo, hdi_lo], [0, yl(2)*0.3], 'r-', 'LineWidth', 1.5);
            plot([hdi_hi, hdi_hi], [0, yl(2)*0.3], 'r-', 'LineWidth', 1.5);

            de_val = [rr.eps_inf_E, rr.eps_inf_L, rr.f_relax/1e9];
            xline(de_val(iP), '--k', 'LineWidth', 1);

            ylabel(param_names{iP}, 'FontName', fontName, 'FontSize', 9);
            if iP == 1
                title(sprintf('%s (R-hat=%.3f)', rr.name, max(rr.rhat)), ...
                    'FontName', fontName, 'FontSize', 10, 'FontWeight', 'bold', ...
                    'Interpreter', 'none');
            end
            set(gca, 'FontName', fontName, 'FontSize', 8, 'Box', 'on');
        end
    end
    sgtitle('Posterior Marginal Densities (95% HDI in red, DE optimum dashed)', ...
        'FontName', fontName, 'FontSize', 14, 'FontWeight', 'bold');
    exportgraphics(fig2, fullfile(outDir, 'Posterior_Density.png'), 'Resolution', 150);
    fprintf('  Saved: Posterior_Density.png\n');

    % ================================================================
    % FIGURE 3: Pairwise Posterior Correlations
    % ================================================================
    fig3 = figure('Color', 'w', 'Position', [50, 50, 1000, 300*n_samples]);

    pair_labels = {'\epsilon_\infty^E', '\epsilon_\infty^L', 'f_{relax} (GHz)'};

    for iF = 1:n_samples
        rr = results(iF);
        posterior = [];
        for c = 1:n_chains
            posterior = [posterior; chains{c}];
        end

        pairs = [1 2; 1 3; 2 3];

        for p = 1:3
            subplot(n_samples, 3, (iF-1)*3 + p);
            p1 = posterior(:, pairs(p,1));
            p2 = posterior(:, pairs(p,2));
            if pairs(p,2) == 3, p2 = p2 / 1e9; end
            if pairs(p,1) == 3, p1 = p1 / 1e9; end

            n_plot = min(5000, size(posterior, 1));
            idx = round(linspace(1, size(posterior, 1), n_plot));
            scatter(p1(idx), p2(idx), 3, [0.3 0.5 0.8], 'filled', ...
                'MarkerFaceAlpha', 0.15);

            xlabel(pair_labels{pairs(p,1)}, 'FontName', fontName, 'FontSize', 10);
            ylabel(pair_labels{pairs(p,2)}, 'FontName', fontName, 'FontSize', 10);
            if p == 1
                title(rr.name, 'FontName', fontName, 'FontSize', 11, ...
                    'FontWeight', 'bold', 'Interpreter', 'none');
            end
            set(gca, 'FontName', fontName, 'FontSize', 9, 'Box', 'on');
        end
    end
    sgtitle('Posterior Pairwise Correlations', ...
        'FontName', fontName, 'FontSize', 14, 'FontWeight', 'bold');
    exportgraphics(fig3, fullfile(outDir, 'Posterior_Correlation.png'), 'Resolution', 150);
    fprintf('  Saved: Posterior_Correlation.png\n');

    % ================================================================
    % FIGURE 4: Posterior Predictive |S21|
    % ================================================================
    fig4 = figure('Color', 'w', 'Position', [50, 50, 1000, 700]);

    for iF = 1:n_samples
        subplot(2, 2, iF);
        hold on;
        rr = results(iF);
        fplot = rr.f_fit / 1e9;

        plot(fplot, 20*log10(abs(rr.S_meas)), '.', ...
            'Color', [0.3 0.3 0.3], 'MarkerSize', 3);

        posterior = [];
        for c = 1:n_chains
            posterior = [posterior; chains{c}];
        end
        n_draws = 100;
        idx_draws = round(linspace(1, size(posterior, 1), n_draws));
        S_pred = zeros(length(fplot), n_draws);
        for d = 1:n_draws
            xd = posterior(idx_draws(d), :);
            Sd = tmm_debye_S21(xd(1), xd(2), xd(3), ...
                delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, ...
                rr.f_fit, c0, Z0);
            S_pred(:, d) = 20*log10(abs(Sd));
        end
        pred_median = median(S_pred, 2);
        pred_lo = quantile(S_pred, 0.025, 2);
        pred_hi = quantile(S_pred, 0.975, 2);

        fill([fplot; flipud(fplot)], [pred_lo; flipud(pred_hi)], ...
            [0.85 0.17 0.15], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
        plot(fplot, pred_median, '-', 'Color', [0.85 0.17 0.15], 'LineWidth', 1.5);

        title(sprintf('%s (RMS=%.2f dB)', rr.name, rr.rms), ...
            'FontName', fontName, 'FontSize', 11, 'FontWeight', 'bold', ...
            'Interpreter', 'none');
        xlabel('Frequency (GHz)', 'FontName', fontName, 'FontSize', 11);
        ylabel('|S_{21}| (dB)', 'FontName', fontName, 'FontSize', 11);
        legend({'Measured', '95% CI', 'Median'}, ...
            'FontName', fontName, 'FontSize', 8, 'Location', 'southwest');
        set(gca, 'FontName', fontName, 'FontSize', 10, 'Box', 'on');
        grid on;
    end
    sgtitle('Posterior Predictive |S_{21}|', ...
        'FontName', fontName, 'FontSize', 14, 'FontWeight', 'bold');
    exportgraphics(fig4, fullfile(outDir, 'Posterior_Predictive.png'), 'Resolution', 150);
    fprintf('  Saved: Posterior_Predictive.png\n');

    % ================================================================
    % FIGURE 5: Parameter Summary Bar Chart
    % ================================================================
    fig5 = figure('Color', 'w', 'Position', [50, 50, 1000, 400]);

    x_labels = cellfun(@(s) strrep(s, '_', '\_'), {results.name}, ...
        'UniformOutput', false);
    param_colors = {[0.29 0.49 0.73], [0.86 0.65 0.23], [0.47 0.67 0.19]};

    for iP = 1:n_params
        subplot(1, 3, iP);
        hold on;

        means = zeros(1, n_samples);
        los = zeros(1, n_samples);
        his = zeros(1, n_samples);
        for iF = 1:n_samples
            rr = results(iF);
            m = rr.post_mean(iP);
            if iP == 3, m = m / 1e9; end
            lo = rr.post_hdi_lo(iP);
            if iP == 3, lo = lo / 1e9; end
            hi = rr.post_hdi_hi(iP);
            if iP == 3, hi = hi / 1e9; end
            means(iF) = m;
            los(iF)  = m - lo;
            his(iF)  = hi - m;
        end

        bar(means, 'FaceColor', param_colors{iP}, 'EdgeColor', 'k', 'LineWidth', 0.5);
        hold on;
        errorbar(1:n_samples, means, los, his, 'k.', 'LineWidth', 1.2, 'CapSize', 8);

        set(gca, 'XTickLabel', x_labels, 'FontName', fontName, 'FontSize', 9, 'Box', 'on');
        xtickangle(30);
        ylabel(param_names{iP}, 'FontName', fontName, 'FontSize', 11);
        title(param_names{iP}, 'FontName', fontName, 'FontSize', 12, 'FontWeight', 'bold');
    end
    sgtitle('Posterior Means with 95% HDI Error Bars', ...
        'FontName', fontName, 'FontSize', 14, 'FontWeight', 'bold');
    exportgraphics(fig5, fullfile(outDir, 'Posterior_Summary.png'), 'Resolution', 150);
    fprintf('  Saved: Posterior_Summary.png\n');

    close all;
end
