# Phase-Enriched Cost + Bayesian MCMC — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add S21 phase to cost function and replace deterministic DE with adaptive MCMC for posterior uncertainty quantification.

**Architecture:** Extract TMM forward model to standalone file, add phase-enriched cost function, implement adaptive Metropolis-Hastings sampler, integrate MCMC after DE in the main pipeline, add posterior visualization. Parallel execution via `parfor` over samples × chains.

**Tech Stack:** MATLAB R2024b, Parallel Computing Toolbox, Statistics and Machine Learning Toolbox, Signal Processing Toolbox

---

## File structure

| File | Action | Purpose |
|------|--------|---------|
| `simulations/tmm_forward.m` | **Create** | Extracted `tmm_debye_S21` + `tmm_chebyshev` (unchanged logic, now standalone) |
| `simulations/cost_tmm3_phase.m` | **Create** | Phase-enriched cost + log-prior + log-likelihood functions |
| `simulations/mcmc_sampler.m` | **Create** | Adaptive Metropolis-Hastings with Gelman-Rubin R-hat |
| `simulations/plot_posterior.m` | **Create** | 5 posterior figures |
| `simulations/run_tmm_inversion.m` | **Modify** | Remove extracted subfunctions, use new files, add MCMC + parfor stage |

### Why extraction is necessary

The TMM forward model (`tmm_debye_S21`, `tmm_chebyshev`) and cost function are currently subfunctions inside `run_tmm_inversion.m`. For MCMC running inside `parfor`, all called functions must be available as separate files on the path. Extracting them also fixes the code duplication noted in CLAUDE.md.

---

### Task 1: Extract TMM forward model to standalone file

**Files:**
- Create: `simulations/tmm_forward.m`
- Modify: `simulations/run_tmm_inversion.m` — remove lines 275-350

- [ ] **Step 1: Create `simulations/tmm_forward.m`**

Copy the two forward-model functions from `run_tmm_inversion.m` lines 275-350 into a new file. The functions are unchanged — only moved.

```matlab
function [A, B, C, D] = tmm_chebyshev(eps_E, eps_L, d_E, d_L, N_pairs, f, c0, Z0)
    % Compute total ABCD matrix using Chebyshev polynomial identity.
    % All inputs are N_freq x 1 vectors (eps_E, eps_L may be complex).
    % d_E, d_L are scalars.
    % Returns A, B, C, D as N_freq x 1 column vectors.

    n_E = sqrt(eps_E);
    n_L = sqrt(eps_L);

    k0 = 2 * pi * f / c0;

    delta_E = k0 .* n_E * d_E;
    delta_L = k0 .* n_L * d_L;

    Z_E = Z0 ./ n_E;
    Z_L = Z0 ./ n_L;

    cE = cos(delta_E); sE = sin(delta_E);
    cL = cos(delta_L); sL = sin(delta_L);

    alpha = cE.*cL - (Z_E./Z_L).*sE.*sL;
    beta  = 1j * (Z_E.*sE.*cL + Z_L.*cE.*sL);
    gamma = 1j * ((1./Z_E).*sE.*cL + (1./Z_L).*cE.*sL);
    delta_term = cE.*cL - (Z_L./Z_E).*sE.*sL;

    s = (alpha + delta_term) / 2;
    theta = acos(s);

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

    A = U_Nm1 .* alpha      - U_Nm2;
    B = U_Nm1 .* beta;
    C = U_Nm1 .* gamma;
    D = U_Nm1 .* delta_term - U_Nm2;
end

function S21 = tmm_debye_S21(eps_inf_E, eps_inf_L, f_relax, ...
        delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, f, c0, Z0)
    % Compute S21 for alternating E-L layered structure.
    % Each layer uses Debye dispersion.
    % Uses Chebyshev acceleration: O(1) per frequency, independent of N_pairs.

    eps_E = eps_inf_E + delta_eps_E ./ (1 + 1j * f / f_relax);
    eps_L = eps_inf_L + delta_eps_L ./ (1 + 1j * f / f_relax);

    [A, B, C, D] = tmm_chebyshev(eps_E, eps_L, d_E, d_L, N_pairs, f, c0, Z0);

    denom = A + B/Z0 + C*Z0 + D;
    S21 = 2 ./ denom;
end
```

- [ ] **Step 2: Remove extracted functions from `run_tmm_inversion.m`**

Delete lines 275-350 in `simulations/run_tmm_inversion.m` (the `tmm_debye_S21` and `tmm_chebyshev` function definitions).

- [ ] **Step 3: Verify extraction — run the existing DE pipeline**

Run in MATLAB:
```matlab
cd('D:/woods/wood-dielectric-inversion');
addpath('simulations');
run_tmm_inversion;
```

Expected: same output as before. The DE should find the same results since `tmm_forward.m` is on the path and contains identical logic.

- [ ] **Step 4: Commit**

```bash
cd D:/woods/wood-dielectric-inversion
git add simulations/tmm_forward.m simulations/run_tmm_inversion.m
git commit -m "refactor: extract TMM forward model to standalone file"
```

---

### Task 2: Create phase-enriched cost function + log-prior + log-likelihood

**Files:**
- Create: `simulations/cost_tmm3_phase.m`

- [ ] **Step 1: Create `simulations/cost_tmm3_phase.m`**

Contains three functions: the cost function (for DE), the log-prior, and the log-likelihood (for MCMC).

```matlab
function J = cost_tmm3_phase(x, f, mag_meas_dB, phase_meas, ...
        delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, c0, Z0)
    % Phase-enriched cost function for DE optimization.
    % x = [eps_inf_E, eps_inf_L, f_relax]
    % J = 0.6 * Huber(|S21| error) + 0.4 * wrapped phase MSE

    S21 = tmm_debye_S21(x(1), x(2), x(3), ...
        delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, f, c0, Z0);

    % Magnitude error (Huber loss)
    mag_err = 20 * log10(max(abs(S21), 1e-12)) - mag_meas_dB;
    L_mag = mag_err.^2;
    mask = abs(mag_err) > 1.0;
    L_mag(mask) = 2.0 * abs(mag_err(mask)) - 1.0;

    % Phase error (wrapped to [-pi, pi])
    phase_model = angle(S21);
    phase_err = atan2(sin(phase_model - phase_meas), ...
                      cos(phase_model - phase_meas));
    L_phase = phase_err.^2;

    J = 0.6 * mean(L_mag) + 0.4 * mean(L_phase);
end

function lp = log_prior_tmm3(x, lb, ub)
    % Log-prior for TMM 3-parameter model.
    % Uniform on eps_inf_E, eps_inf_L within [lb(1:2), ub(1:2)].
    % Log-uniform on f_relax within [lb(3), ub(3)]: p(x) ~ 1/x.

    if any(x < lb) || any(x > ub)
        lp = -inf;
        return;
    end
    lp = -log(x(3));  % log-uniform prior on f_relax
end

function ll = log_like_tmm3(x, f, mag_meas_dB, phase_meas, ...
        delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, c0, Z0)
    % Log-likelihood for MCMC.
    % ll = -J, where J is the phase-enriched cost.
    % Constants (sigma, 2*pi) cancel in M-H acceptance ratio.

    J = cost_tmm3_phase(x, f, mag_meas_dB, phase_meas, ...
        delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, c0, Z0);
    ll = -J;
end
```

- [ ] **Step 2: Verify cost function compiles and runs**

Run in MATLAB:
```matlab
cd('D:/woods/wood-dielectric-inversion');
addpath('simulations');
c0 = 299792458; Z0 = 377;
f = linspace(18e9, 44e9, 100)';
mag_dB = -30 * ones(100, 1);
phase = zeros(100, 1);
x = [3.0, 4.0, 2e9];
J = cost_tmm3_phase(x, f, mag_dB, phase, 3.0, 1.5, 3e-3, 2e-3, 24, c0, Z0);
disp(J);
```
Expected: a finite scalar (not NaN, not Inf).

- [ ] **Step 3: Commit**

```bash
git add simulations/cost_tmm3_phase.m
git commit -m "feat: add phase-enriched cost function with log-prior and log-likelihood"
```

---

### Task 3: Create adaptive Metropolis-Hastings MCMC sampler

**Files:**
- Create: `simulations/mcmc_sampler.m`

- [ ] **Step 1: Create `simulations/mcmc_sampler.m`**

```matlab
function [chains, diagnostics] = mcmc_sampler(log_prior_fn, log_like_fn, ...
        init_x, lb, ub, n_burnin, n_samples, adapt_interval)
    % Adaptive Metropolis-Hastings MCMC (Haario et al., 2001).
    %
    % Inputs:
    %   log_prior_fn   — function handle: lp = f(x, lb, ub)
    %   log_like_fn    — function handle: ll = f(x)
    %   init_x         — n_chains x D starting points
    %   lb, ub         — 1 x D bounds
    %   n_burnin       — burn-in iterations per chain
    %   n_samples      — post-burn-in samples per chain
    %   adapt_interval — update proposal covariance every N steps
    %
    % Outputs:
    %   chains      — cell array {n_chains}(n_samples, D)
    %   diagnostics — struct with .rhat (1xD) and .acceptance (n_chainsx1)

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
                % Scale factor from Gelman et al. (1997): (2.38^2)/D
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
```

- [ ] **Step 2: Quick unit test with a toy Gaussian posterior**

Run in MATLAB:
```matlab
cd('D:/woods/wood-dielectric-inversion');
addpath('simulations');

% True: 2-D Gaussian at (3, 5)
log_prior = @(x, lb, ub) 0;
log_like  = @(x) -0.5 * sum((x - [3, 5]).^2) / 2^2;

init = [2, 4; 4, 6; 2.5, 5.5; 3.5, 4.5];
lb = [-10, -10];
ub = [10, 10];

[chains, diag] = mcmc_sampler(log_prior, log_like, init, lb, ub, 2000, 5000, 200);
post = cell2mat(cellfun(@(c) c, chains(:), 'UniformOutput', false));

fprintf('True mean: [3, 5]\n');
fprintf('Posterior mean: [%.2f, %.2f]\n', mean(post));
fprintf('R-hat: [%.4f, %.4f]\n', diag.rhat);
fprintf('Acceptance rate: [%.2f, %.2f, %.2f, %.2f]\n', diag.acceptance);
```

Expected: posterior mean close to (3, 5), R-hat < 1.1, acceptance ~15-40%.

- [ ] **Step 3: Commit**

```bash
git add simulations/mcmc_sampler.m
git commit -m "feat: add adaptive Metropolis-Hastings MCMC sampler"
```

---

### Task 4: Integrate phase cost + MCMC into main pipeline

**Files:**
- Modify: `simulations/run_tmm_inversion.m`

- [ ] **Step 1: Modify the DE optimizer subfunction to use phase-enriched cost**

In `simulations/run_tmm_inversion.m`, replace the `cost_tmm3` call in `tmm_debye_invert` (around line 368) with `cost_tmm3_phase`. Also save phase data for the cost function.

Replace the `tmm_debye_invert` function (lines 352-412) with a version that accepts and passes phase data:

```matlab
function [best_x, loss_curve] = tmm_debye_invert(S_meas, f, lb, ub, ...
        delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, ...
        popSize, maxIter, c0, Z0)
    % Differential Evolution for 3-parameter TMM+Debye model.
    % Uses phase-enriched cost function.

    D = 3;
    mag_meas_dB = 20 * log10(max(abs(S_meas), 1e-12));
    phase_meas  = angle(S_meas);

    % Latin Hypercube initialization
    X = repmat(lb, popSize, 1) + lhsdesign(popSize, D) .* repmat(ub - lb, popSize, 1);

    cost = zeros(popSize, 1);
    for i = 1:popSize
        cost(i) = cost_tmm3_phase(X(i,:), f, mag_meas_dB, phase_meas, ...
            delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, c0, Z0);
    end

    [best_cost, idx] = min(cost);
    best_x = X(idx, :);
    loss_curve = zeros(maxIter, 1);

    F  = 0.5;
    CR = 0.7;

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

            cost_u = cost_tmm3_phase(u, f, mag_meas_dB, phase_meas, ...
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
```

Also remove the old `cost_tmm3` subfunction (lines 414-429).

- [ ] **Step 2: Add MCMC stage after DE in the main loop**

In `run_tmm_inversion.m`, after the DE call (~line 129), add MCMC. Replace the block from the DE call to the result storage with:

```matlab
        % ---- DE Optimization (3-param TMM + Debye, phase-enriched) ----
        rng(iF * 100);
        [p_de, loss_curve] = tmm_debye_invert(...
            S_fit_sm, f_fit, lb3, ub3, ...
            delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, ...
            popSize, maxIter, c0, Z0);

        % ---- MCMC: 4-chain adaptive Metropolis-Hastings ----
        n_chains  = 4;
        n_burnin  = 5000;
        n_samples = 20000;
        adapt_int = 500;

        % Initialize 4 chains from DE best-fit with perturbation
        rng_state = rng;
        init_x = repmat(p_de, n_chains, 1) .* ...
            (1 + 0.02 * randn(n_chains, 3));  % 2% perturbation
        init_x = max(lb3, min(ub3, init_x));
        rng(rng_state);  % restore for reproducibility

        % Prepare function handles for MCMC
        mag_meas_dB = 20 * log10(max(abs(S_fit_sm), 1e-12));
        phase_meas  = angle(S_fit_sm);

        % Build anonymous functions binding data
        lp_fn = @(x, l, u) log_prior_tmm3(x, l, u);
        ll_fn = @(x) log_like_tmm3(x, f_fit, mag_meas_dB, phase_meas, ...
            delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, c0, Z0);

        [mcmc_chains, mcmc_diag] = mcmc_sampler(lp_fn, ll_fn, ...
            init_x, lb3, ub3, n_burnin, n_samples, adapt_int);

        % Pool all posterior samples for summary stats
        posterior = [];
        for c = 1:n_chains
            posterior = [posterior; mcmc_chains{c}];
        end

        % ---- Store ----
        eps_inf_E = p_de(1);
        eps_inf_L = p_de(2);
        f_relax   = p_de(3);

        S21_model = tmm_debye_S21(eps_inf_E, eps_inf_L, f_relax, ...
            delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, f_fit, c0, Z0);
        rms_val = sqrt(mean((20*log10(abs(S21_model)) - 20*log10(abs(S_fit_sm))).^2));

        results(iF).name       = baseName;
        results(iF).eps_inf_E  = eps_inf_E;
        results(iF).eps_inf_L  = eps_inf_L;
        results(iF).f_relax    = f_relax;
        results(iF).rms        = rms_val;
        results(iF).f_fit      = f_fit;
        results(iF).S_meas     = S_fit_sm;
        results(iF).S_model    = S21_model;
        results(iF).loss_curve = loss_curve;

        % MCMC results
        results(iF).post_mean  = mean(posterior, 1);
        hdi_lo = quantile(posterior, 0.025, 1);
        hdi_hi = quantile(posterior, 0.975, 1);
        results(iF).post_hdi_lo = hdi_lo;
        results(iF).post_hdi_hi = hdi_hi;
        results(iF).rhat       = mcmc_diag.rhat;
        results(iF).accept_rate = mcmc_diag.acceptance;
        results(iF).mcmc_chains = mcmc_chains;
```

- [ ] **Step 3: Update console output to include MCMC summary**

After the main loop, add posterior summary printing. Insert before the Figure 1 code:

```matlab
    % =================================================================
    % MCMC Posterior Summary
    % =================================================================
    fprintf('\n=============== MCMC Posterior Summary ===============\n');
    fprintf('%-22s %-28s %-28s %-28s %6s\n', ...
        'Sample', 'eps_inf_E (95%% HDI)', 'eps_inf_L (95%% HDI)', ...
        'f_relax GHz (95%% HDI)', 'R-hat');
    fprintf('%s\n', repmat('-', 1, 115));
    for iF = 1:length(results)
        rr = results(iF);
        fprintf('%-22s %6.2f [%5.2f, %5.2f]  %6.2f [%5.2f, %5.2f]  %5.1f [%4.1f, %4.1f]    %5.3f\n', ...
            rr.name, ...
            rr.post_mean(1), rr.post_hdi_lo(1), rr.post_hdi_hi(1), ...
            rr.post_mean(2), rr.post_hdi_lo(2), rr.post_hdi_hi(2), ...
            rr.post_mean(3)/1e9, rr.post_hdi_lo(3)/1e9, rr.post_hdi_hi(3)/1e9, ...
            max(rr.rhat));
        fprintf('    Acceptance: [%.2f %.2f %.2f %.2f]\n', rr.accept_rate);
    end

    fprintf('\n--- Texture Anisotropy (Posterior Means) ---\n');
    for iPair = 1:2
        vert_name = pairs_cfg{iPair}{1};
        para_name = pairs_cfg{iPair}{2};
        vert_mean = []; para_mean = [];
        for iF = 1:length(results)
            if strcmp(results(iF).name, vert_name)
                vert_mean = results(iF).post_mean;
            end
            if strcmp(results(iF).name, para_name)
                para_mean = results(iF).post_mean;
            end
        end
        if ~isempty(vert_mean) && ~isempty(para_mean)
            fprintf('%s:\n', pair_titles{iPair});
            fprintf('  eps_inf_E: %.2f vs %.2f (delta=%.2f)\n', ...
                vert_mean(1), para_mean(1), vert_mean(1)-para_mean(1));
            fprintf('  eps_inf_L: %.2f vs %.2f (delta=%.2f)\n', ...
                vert_mean(2), para_mean(2), vert_mean(2)-para_mean(2));
            fprintf('  f_relax:   %.1f vs %.1f GHz (delta=%.1f)\n', ...
                vert_mean(3)/1e9, para_mean(3)/1e9, (vert_mean(3)-para_mean(3))/1e9);
        end
    end
```

- [ ] **Step 4: Add plot_posterior call after existing figures**

Append before the `save` call at the end of the script:

```matlab
    % =================================================================
    % Posterior Visualization
    % =================================================================
    plot_posterior(results, outDir, fontName, ...
        delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, c0, Z0);
```

- [ ] **Step 5: Update the save call**

Change the final `save` call to include MCMC data:

```matlab
    save(fullfile(outDir, 'MCMC_Inversion_Results.mat'), 'results', ...
        'd_total', 'd_E_mm', 'd_L_mm', 'N_pairs', 'delta_eps_E', 'delta_eps_L', ...
        'n_burnin', 'n_samples', 'n_chains');
    fprintf('\nResults saved to: %s\n', fullfile(outDir, 'MCMC_Inversion_Results.mat'));
```

- [ ] **Step 6: Verify DE stage still works**

Run in MATLAB:
```matlab
cd('D:/woods/wood-dielectric-inversion');
addpath('simulations');
run_tmm_inversion;
```

Expected: DE optimization completes without errors, parameter values may shift slightly from before (phase now contributes to cost). Console shows both DE and MCMC output.

- [ ] **Step 7: Commit**

```bash
git add simulations/run_tmm_inversion.m
git commit -m "feat: integrate phase-enriched cost and MCMC into TMM inversion pipeline"
```

---

### Task 5: Create posterior visualization

**Files:**
- Create: `simulations/plot_posterior.m`

- [ ] **Step 1: Create `simulations/plot_posterior.m`**

```matlab
function plot_posterior(results, outDir, fontName, ...
        delta_eps_E, delta_eps_L, d_E, d_L, N_pairs, c0, Z0)
    % Generate 5 posterior diagnostic figures.

    n_samples = length(results);
    n_params = 3;
    param_names = {'\epsilon_\infty^E', '\epsilon_\infty^L', 'f_{relax} (GHz)'};
    sample_names = {results.name};
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
                if iP == 3, y = y / 1e9; end  % f_relax in GHz
                plot(y, 'Color', [chain_colors(c,:) 0.4], 'LineWidth', 0.3);
            end
            ylabel(param_names{iP}, 'FontName', fontName, 'FontSize', 9);
            if iP == 1
                title(sample_names{iF}, 'FontName', fontName, 'FontSize', 10, ...
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

            % Kernel density estimate
            [f_est, xi] = ksdensity(vals);
            fill([xi, fliplr(xi)], [f_est, zeros(size(f_est))], ...
                [0.7 0.85 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.6);

            % 95% HDI
            hdi_lo = quantile(vals, 0.025);
            hdi_hi = quantile(vals, 0.975);
            yl = ylim;
            plot([hdi_lo, hdi_lo], [0, yl(2)*0.3], 'r-', 'LineWidth', 1.5);
            plot([hdi_hi, hdi_hi], [0, yl(2)*0.3], 'r-', 'LineWidth', 1.5);

            % DE best-fit marker
            de_val = [rr.eps_inf_E, rr.eps_inf_L, rr.f_relax/1e9];
            xline(de_val(iP), '--k', 'LineWidth', 1);

            ylabel(param_names{iP}, 'FontName', fontName, 'FontSize', 9);
            if iP == 1
                title(sprintf('%s (R-hat=%.3f)', sample_names{iF}, max(rr.rhat)), ...
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

    for iF = 1:n_samples
        rr = results(iF);
        posterior = [];
        for c = 1:n_chains
            posterior = [posterior; chains{c}];
        end

        pair_labels = {'\epsilon_\infty^E', '\epsilon_\infty^L', 'f_{relax} (GHz)'};
        pairs = [1 2; 1 3; 2 3];

        for p = 1:3
            subplot(n_samples, 3, (iF-1)*3 + p);
            p1 = posterior(:, pairs(p,1));
            p2 = posterior(:, pairs(p,2));
            if pairs(p,2) == 3, p2 = p2 / 1e9; end
            if pairs(p,1) == 3, p1 = p1 / 1e9; end

            % Thin for plot performance
            n_plot = min(5000, size(posterior, 1));
            idx = round(linspace(1, size(posterior, 1), n_plot));
            scatter(p1(idx), p2(idx), 3, [0.3 0.5 0.8], 'filled', ...
                'MarkerFaceAlpha', 0.15);

            xlabel(pair_labels{pairs(p,1)}, 'FontName', fontName, 'FontSize', 10);
            ylabel(pair_labels{pairs(p,2)}, 'FontName', fontName, 'FontSize', 10);
            if p == 1
                title(sample_names{iF}, 'FontName', fontName, 'FontSize', 11, ...
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

        % Measured data
        plot(fplot, 20*log10(abs(rr.S_meas)), '.', ...
            'Color', [0.3 0.3 0.3], 'MarkerSize', 3);

        % Posterior predictive: 100 draws
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
                fplot*1e9, c0, Z0);
            S_pred(:, d) = 20*log10(abs(Sd));
        end
        pred_median = median(S_pred, 2);
        pred_lo = quantile(S_pred, 0.025, 2);
        pred_hi = quantile(S_pred, 0.975, 2);

        % 95% credible band
        fill([fplot; flipud(fplot)], [pred_lo; flipud(pred_hi)], ...
            [0.85 0.17 0.15], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
        plot(fplot, pred_median, '-', 'Color', [0.85 0.17 0.15], 'LineWidth', 1.5);

        rms_val = rr.rms;
        title(sprintf('%s (RMS=%.2f dB)', rr.name, rms_val), ...
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

    x_labels = {results.name};
    x_labels = cellfun(@(s) strrep(s, '_', '\_'), x_labels, 'UniformOutput', false);

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

        % Group bars: vertical first, then parallel
        b = bar(means, 'FaceColor', param_colors{iP}, 'EdgeColor', 'k', 'LineWidth', 0.5);
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
```

- [ ] **Step 2: Test the visualization functions with existing data**

If there are already results from a previous DE run, skip. Otherwise, run the full pipeline after Task 4 is complete.

- [ ] **Step 3: Commit**

```bash
git add simulations/plot_posterior.m
git commit -m "feat: add posterior visualization (5 diagnostic figures)"
```

---

### Task 6: Parallel execution — add parfor

**Files:**
- Modify: `simulations/run_tmm_inversion.m`

- [ ] **Step 1: Wrap the sample loop with parfor and add pool setup**

Add at the top of `run_tmm_inversion` (before the sample loop), replace the `for iF = 1:length(matFiles)` with `parfor`:

```matlab
    % ====== Start Parallel Pool ======
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', min(16, feature('numCores')));
    end

    % ... [frequency setup, constants remain unchanged] ...

    % ====== Process Each Sample (PARALLEL) ======
    parfor iF = 1:length(matFiles)
        % ... [full sample processing block from Task 4] ...
    end
```

Note: `parfor` requires that all variables used inside the loop are available. The `results` struct must be constructed to support `parfor` slicing — use indexed assignment on the loop variable.

- [ ] **Step 2: Ensure parfor compatibility**

Move the `results` struct initialization to pre-allocate before `parfor`:

```matlab
    results = struct();
    % Pre-allocate for parfor compatibility
    results_arr = cell(length(matFiles), 1);
```

Replace `results(iF).field = value` with:
```matlab
    tmp = struct();
    tmp.name = baseName;
    tmp.eps_inf_E = eps_inf_E;
    % ... all fields ...
    results_arr{iF} = tmp;
```

And after the parfor loop, convert back:
```matlab
    results = [results_arr{:}];
```

- [ ] **Step 3: Full pipeline test with parallel execution**

Run in MATLAB:
```matlab
cd('D:/woods/wood-dielectric-inversion');
addpath('simulations');
run_tmm_inversion;
```

Expected: all 4 samples processed in parallel, MCMC chains run, 5 figures saved to `simulations/results/`, .mat file saved. Runtime < 20 minutes.

- [ ] **Step 4: Verify output files exist**

```matlab
dir('simulations/results/Posterior_*.png')
dir('simulations/results/MCMC_Inversion_Results.mat')
```

Expected: 5 PNG files + 1 MAT file.

- [ ] **Step 5: Commit**

```bash
git add simulations/run_tmm_inversion.m
git commit -m "feat: add parfor parallel execution for MCMC chains"
```

---

### Task 7: End-to-end validation and cleanup

**Files:**
- No new files. Run and verify.

- [ ] **Step 1: Full clean run**

```matlab
cd('D:/woods/wood-dielectric-inversion');
addpath('simulations');
delete(gcp('nocreate'));
run_tmm_inversion;
```

- [ ] **Step 2: Check convergence diagnostics**

Verify in console output:
- All R-hat values < 1.1
- Acceptance rates between 15% and 50%
- Trace plots show good mixing (no obvious trends)

- [ ] **Step 3: Check texture anisotropy results**

Verify in console output:
- Parallel-texture eps systematically higher than perpendicular
- f_relax in 1-4 GHz range (physically plausible for bound-water Debye)

- [ ] **Step 4: Commit final state**

```bash
git add -A
git status
git commit -m "chore: final validation of phase+MCMC pipeline"
```
