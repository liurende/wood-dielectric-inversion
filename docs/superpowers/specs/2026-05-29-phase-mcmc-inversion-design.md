# Phase-Enriched Cost Function + Bayesian MCMC Inversion

## Summary

Two-part improvement to the TMM + Debye wood dielectric inversion:

1. Add S21 phase to the cost function (currently magnitude-only) — free information gain
2. Replace deterministic DE point estimation with Bayesian MCMC for posterior uncertainty quantification

Hardware: 14-core / 20-thread CPU, MATLAB Parallel Computing Toolbox R2024b available.

## Architecture

```
CSV -> preprocess (unchanged) -> DE global search (new cost: mag + phase)
                                       |
                                       v
                                DE best solution -> MCMC initial point
                                       |
                                       v
                          +-- Chain 1 (worker 1) --+
     parfor over samples: +-- Chain 2 (worker 2) --+--> merge posteriors
                          +-- Chain 3 (worker 3) --+         |
                          +-- Chain 4 (worker 4) --+         v
                                              post-processing: trace plots,
                                              posterior densities, correlation
                                              pairs, predictive bands, summary
```

### Files

| File | Action | Purpose |
|------|--------|---------|
| `simulations/run_tmm_inversion.m` | Modify | Add phase to cost, add MCMC stage after DE |
| `simulations/mcmc_sampler.m` | **New** | Adaptive Metropolis-Hastings with convergence diagnostics |
| `simulations/plot_posterior.m` | **New** | Posterior visualization (5 figures) |
| `simulations/results/` | Existing | All output figures and .mat files saved here |

Preprocessing, TMM forward model (`tmm_debye_S21`, `tmm_chebyshev`), and `load_measurement.m` are unchanged.

## Cost function

Current (magnitude-only):
```
J_mag = mean(Huber(|S21_model|_dB - |S21_meas|_dB))
```

New (magnitude + phase):
```
J = 0.6 * J_mag + 0.4 * J_phase
```

Phase error handles 2π wrapping:
```
phase_err = angle(S21_model) - angle(S21_meas)
phase_err = atan2(sin(phase_err), cos(phase_err))  % wrap to [-π, π]
J_phase = mean(phase_err.^2)
```

Weight rationale: magnitude has higher SNR, but phase is more sensitive to thickness and f_relax — it helps break the remaining parameter coupling.

## MCMC setup

- **Algorithm**: Adaptive Metropolis-Hastings (Haario et al., 2001)
- **Proposal**: Multivariate Gaussian, covariance adapted from DE population covariance scaled by (2.38^2)/D
- **Chains per sample**: 4
- **Burn-in**: 5000 samples per chain
- **Samples**: 20000 per chain (80,000 total per sample, 20,000 post-burn-in)
- **Convergence**: Gelman-Rubin R-hat < 1.1 across all 4 chains
- **Parallelization**: `parfor` over samples (4) × chains per sample (4) — potentially 16 workers

### Priors

| Parameter | Prior | Range |
|-----------|-------|-------|
| eps_inf_E | Uniform | [1.5, 10.0] |
| eps_inf_L | Uniform | [1.5, 10.0] |
| f_relax | Log-uniform | [0.5e9, 60e9] |

Log-uniform on f_relax: relaxes over orders of magnitude, physically appropriate for relaxation frequencies spanning 1–60 GHz.

### Adaptive proposal

After an initial non-adaptive phase (first 1000 samples), the proposal covariance is updated every 500 samples using the empirical covariance of the chain so far, scaled by (2.38^2)/D where D=3 is the parameter dimension. This achieves near-optimal acceptance rate (~23% for D=3).

## Post-processing outputs

All saved to `simulations/results/`:

1. **`Posterior_Trace.png`** — 4×3 grid of trace plots (4 samples, 3 parameters), showing all 4 chains overlaid per panel. Validates mixing and convergence visually.

2. **`Posterior_Density.png`** — 4×3 grid of posterior marginal densities (KDE) with 95% HDI shaded. DE best-fit value marked as vertical dashed line for comparison.

3. **`Posterior_Correlation.png`** — Pairwise scatter plots of posterior samples (eps_inf_E vs eps_inf_L, eps_inf_E vs f_relax, eps_inf_L vs f_relax), one set per sample. Reveals parameter coupling structure.

4. **`Posterior_Predictive.png`** — 2×2 grid (one per sample). |S21| measurement points overlaid with posterior predictive median and 95% credible band (100 posterior draws). RMS labeled.

5. **`Posterior_Summary.png`** — Bar chart comparing posterior means with 95% HDI error bars across all 4 samples, grouped by parameter. Vertical/parallel texture pairs side by side.

### Console output

Summary table printed to console:
```
Sample              eps_inf_E (95% HDI)    eps_inf_L (95% HDI)    f_relax GHz (95% HDI)    R-hat
LargeWood_Vertical  1.65 [1.52, 1.79]      1.80 [1.68, 1.95]      1.0 [0.8, 1.3]          1.02
...
```

### Saved .mat file

`simulations/results/MCMC_Inversion_Results.mat` containing:
- `chains`: cell array of N_samples×3 posterior sample arrays per chain
- `results`: summary struct with posterior means, HDIs, R-hat per sample
- `mcmc_params`: burn-in, sample count, proposal settings

## Estimated runtime

With 14-core parallel pool:
- DE stage: ~2 min (unchanged, minor overhead from phase calculation)
- MCMC stage: ~10-15 min (16 chains × 25000 evaluations, each ~300 freq points)
- Plotting: ~30 sec
- **Total: ~12-18 minutes**

## Scope note

- S11 joint inversion is NOT included (requires lens antenna hardware not yet available)
- Double Debye model is NOT included (deferred until posterior analysis reveals whether single Debye is the limiting factor)
- The old 4-parameter power-law pipeline (`run_wood_inversion.m` in root) is not modified
