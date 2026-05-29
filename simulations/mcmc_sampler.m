function [chains, diagnostics] = mcmc_sampler(log_prior_fn, log_like_fn, ...
        init_x, lb, ub, n_burnin, n_samples, adapt_interval)
    % Adaptive Metropolis-Hastings MCMC (Haario et al., 2001).
    %
    % Inputs:
    %   log_prior_fn   -- function handle: lp = f(x, lb, ub)
    %   log_like_fn    -- function handle: ll = f(x)
    %   init_x         -- n_chains x D starting points
    %   lb, ub         -- 1 x D bounds
    %   n_burnin       -- burn-in iterations per chain
    %   n_samples      -- post-burn-in samples per chain
    %   adapt_interval -- update proposal covariance every N steps
    %
    % Outputs:
    %   chains      -- cell array {n_chains}(n_samples, D)
    %   diagnostics -- struct with .rhat (1xD), .acceptance_burnin, .acceptance_sampling

    [n_chains, D] = size(init_x);
    assert(n_chains >= 2, 'mcmc_sampler requires at least 2 chains for R-hat diagnostics.');
    total_iter = n_burnin + n_samples;
    chains = cell(n_chains, 1);
    acceptance_burnin  = zeros(n_chains, 1);
    acceptance_sampling = zeros(n_chains, 1);

    for c = 1:n_chains
        x_curr = init_x(c, :);
        lp_curr = log_prior_fn(x_curr, lb, ub);

        % If starting point is invalid, find a valid one
        if ~isfinite(lp_curr)
            for tries = 1:2000
                x_curr = lb + rand(1, D) .* (ub - lb);
                lp_curr = log_prior_fn(x_curr, lb, ub);
                if isfinite(lp_curr), break; end
            end
            if ~isfinite(lp_curr)
                error('Could not find valid starting point for chain %d', c);
            end
        end

        ll_curr = log_like_fn(x_curr);
        chain = zeros(total_iter, D);
        chain(1, :) = x_curr;

        % Initial proposal: diagonal covariance (small steps)
        prop_cov = diag(((ub - lb) * 0.01).^2);
        n_accept_burnin  = 0;
        n_accept_sampling = 0;

        for t = 2:total_iter
            % Propose from multivariate Gaussian
            x_prop = mvnrnd(x_curr, prop_cov);

            % Metropolis-Hastings (out-of-bounds proposals rejected
            % naturally by log-prior returning -inf — preserves symmetry)
            lp_prop = log_prior_fn(x_prop, lb, ub);
            if isfinite(lp_prop)
                ll_prop = log_like_fn(x_prop);
                log_alpha = lp_prop + ll_prop - lp_curr - ll_curr;
                if log_alpha > 0 || log(rand) < log_alpha
                    x_curr = x_prop;
                    lp_curr = lp_prop;
                    ll_curr = ll_prop;
                    if t <= n_burnin
                        n_accept_burnin = n_accept_burnin + 1;
                    else
                        n_accept_sampling = n_accept_sampling + 1;
                    end
                end
            end

            chain(t, :) = x_curr;

            % Adapt proposal covariance during burn-in
            if t <= n_burnin && t > 2*D && mod(t, adapt_interval) == 0
                idx = max(1, t - 500):t;
                chain_sofar = chain(idx, :);
                C = cov(chain_sofar);
                prop_cov = (2.38^2 / D) * C + 1e-8 * eye(D);
            end
        end

        % Discard burn-in
        chains{c} = chain((n_burnin + 1):end, :);
        acceptance_burnin(c)  = n_accept_burnin / (n_burnin - 1);
        acceptance_sampling(c) = n_accept_sampling / n_samples;
    end

    % Split-R-hat diagnostic (Vehtari et al., 2021)
    diagnostics = compute_split_rhat(chains);
    diagnostics.acceptance_burnin  = acceptance_burnin;
    diagnostics.acceptance_sampling = acceptance_sampling;
end

function diagnostics = compute_split_rhat(chains)
    % Split-R-hat (Vehtari et al., 2021).
    % Splits each chain into two halves before computing R-hat.
    % This detects within-chain non-stationarity that original R-hat can miss.

    n_chains = length(chains);
    n_samples = size(chains{1}, 1);
    D = size(chains{1}, 2);

    % Split each chain into two halves
    n_split = 2 * n_chains;
    half = floor(n_samples / 2);
    split_samples = cell(n_split, 1);
    for c = 1:n_chains
        split_samples{2*c-1} = chains{c}(1:half, :);
        split_samples{2*c}   = chains{c}(end-half+1:end, :);
    end

    rhat = zeros(1, D);
    for d = 1:D
        chain_means = zeros(n_split, 1);
        chain_vars  = zeros(n_split, 1);
        for s = 1:n_split
            chain_means(s) = mean(split_samples{s}(:, d));
            chain_vars(s)  = var(split_samples{s}(:, d));
        end
        B = half * var(chain_means);
        W = mean(chain_vars);
        V = (half - 1)/half * W + B/half;
        rhat(d) = sqrt(V / W);
    end
    diagnostics = struct('rhat', rhat);
end
