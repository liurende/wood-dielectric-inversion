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
    %   diagnostics -- struct with .rhat (1xD) and .acceptance (n_chainsx1)

    [n_chains, D] = size(init_x);
    total_iter = n_burnin + n_samples;
    chains = cell(n_chains, 1);
    acceptance = zeros(n_chains, 1);

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
        n_accept = 0;

        for t = 2:total_iter
            % Propose from multivariate Gaussian
            x_prop = mvnrnd(x_curr, prop_cov);

            % Bounce off boundaries
            below = x_prop < lb;
            above = x_prop > ub;
            x_prop(below) = 2 * lb(below) - x_prop(below);
            x_prop(above) = 2 * ub(above) - x_prop(above);
            x_prop = min(ub, max(lb, x_prop));

            % Metropolis-Hastings
            lp_prop = log_prior_fn(x_prop, lb, ub);
            if isfinite(lp_prop)
                ll_prop = log_like_fn(x_prop);
                log_alpha = lp_prop + ll_prop - lp_curr - ll_curr;
                if log_alpha > 0 || log(rand) < log_alpha
                    x_curr = x_prop;
                    lp_curr = lp_prop;
                    ll_curr = ll_prop;
                    n_accept = n_accept + 1;
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
        acceptance(c) = n_accept / total_iter;
    end

    % Gelman-Rubin R-hat diagnostic
    diagnostics = compute_rhat(chains);
    diagnostics.acceptance = acceptance;
end

function diagnostics = compute_rhat(chains)
    % Gelman-Rubin potential scale reduction factor.
    n_chains = length(chains);
    n_samples = size(chains{1}, 1);
    D = size(chains{1}, 2);

    rhat = zeros(1, D);
    for d = 1:D
        chain_means = zeros(n_chains, 1);
        chain_vars  = zeros(n_chains, 1);
        for c = 1:n_chains
            chain_means(c) = mean(chains{c}(:, d));
            chain_vars(c)  = var(chains{c}(:, d));
        end
        B = n_samples * var(chain_means);   % between-chain
        W = mean(chain_vars);               % within-chain
        V = (n_samples - 1)/n_samples * W + B/n_samples;
        rhat(d) = sqrt(V / W);
    end
    diagnostics = struct('rhat', rhat);
end
