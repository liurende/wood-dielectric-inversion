# Wood Dielectric Inversion

Microwave dielectric characterization of wood — inverting complex permittivity from Keysight VNA S₂₁ transmission measurements (18–44 GHz) to study anisotropy between parallel and perpendicular grain orientations.

## Overview

Wood is a natural anisotropic dielectric. When an EM wave propagates through it, the attenuation and phase delay depend on grain orientation relative to the E-field. This project applies global optimization and Bayesian MCMC to extract the complex permittivity tensor from free-space transmission measurements.

Two forward-model pipelines coexist:

| Pipeline | Forward model | Dispersion | Parameters | Status |
|----------|-------------|------------|-----------|--------|
| 1. Power-law | Fabry-Perot | tan δ(f) = tan δ_ref · (f/30GHz)^n | 4 | Baseline — proven underdetermined |
| 2. TMM + Debye | Chebyshev-accelerated transfer matrix | ε* = ε_∞ + Δε/(1+jf/f_relax) | 3 | **Preferred** — well-posed |

Pipeline 2 (TMM + Debye) is the current state of the art. It models wood as alternating earlywood/latewood layers with Debye relaxation, uses physically motivated parameter reductions to achieve well-posedness, and provides full Bayesian posterior distributions via adaptive MCMC.

## Directory Structure

```
.
├── run_wood_inversion.m          # Pipeline 1: 5-algorithm batch inversion (power-law model)
├── plot_wood_comparison.m        # Pipeline 1: 4-figure comparison from saved .mat files
├── n_scan_diagnosis.m            # Pipeline 1: diagnose loss exponent n with fmincon
├── load_measurement.m            # Keysight VNA CSV parser
├── lzmwoods/                     # Raw measurement data
│   ├── air1.csv, air2.csv, air3.csv          # Air reference (through-calibration)
│   ├── 0_垂直.csv, 0_平行.csv                 # Large wood block (vertical/parallel grain)
│   ├── 1_垂直.csv, 1_平行.csv                 # Small wood block (~120 mm)
│   └── results_Wood/                         # Pipeline 1 inversion outputs
├── simulations/                  # Pipeline 2: TMM + Debye + MCMC
│   ├── run_tmm_inversion.m       # Main inversion — loads data, runs DE+MCMC, generates figures
│   ├── tmm_debye_S21.m           # TMM forward model (Chebyshev-accelerated, O(1) per freq)
│   ├── tmm_derivation.md         # Full mathematical derivation (sections 1–9)
│   ├── cost_tmm3_phase.m         # Phase-enriched cost function
│   ├── log_prior_tmm3.m          # Uniform log-prior with bounds
│   ├── log_like_tmm3.m           # Laplace log-likelihood (robust to outliers)
│   ├── mcmc_sampler.m            # Adaptive Metropolis-Hastings MCMC (Haario et al., 2001)
│   ├── plot_posterior.m          # 5-figure posterior diagnostic suite
│   ├── exp1_identifiability.m    # Parameter identifiability — cost landscape analysis
│   ├── exp2_recovery_test.m      # Monte Carlo recovery test under noise
│   ├── exp3_model_mismatch.m     # Debye→power-law model mismatch diagnosis
│   ├── exp4_debye_model.m        # Debye model forward validation
│   ├── exp5_fixed_d.m            # Constrained inversion with known thickness
│   ├── exp6_tmm_layered.m        # Layered vs uniform TMM comparison + contrast sweep
│   ├── group_meeting_figures.m   # 6 presentation-quality figures
│   ├── synthesis.m               # Loads experiment outputs, prints summary
│   └── results/                  # Simulation experiment outputs + figures
└── docs/                         # Design documents (specs, plans)
```

## Physical Models

### Pipeline 1 — Power-Law Loss Tangent + Fabry-Perot

Uses a uniform-medium Fabry-Perot transmission formula with a power-law frequency-dependent loss tangent:

$$\varepsilon^*(f) = \varepsilon_r \cdot \big(1 - j \cdot \tan\delta(f)\big)$$

$$\tan\delta(f) = \tan\delta_{ref} \cdot \left(\frac{f}{30 \text{ GHz}}\right)^n$$

Optimization parameters: `[ε_r, log₁₀(tanδ_ref), n, d]` (4-dimensional).

**Known issue:** (tanδ, n) and (ε_r, d) pairs are fundamentally coupled. The cost landscape has a long narrow valley — these parameters cannot be reliably separated from |S₂₁| magnitude alone. This pipeline is retained for comparison but is **not recommended** as the primary workflow.

### Pipeline 2 — TMM + Debye Relaxation (Preferred)

Models wood as alternating earlywood (E) / latewood (L) layers, each with Debye dispersion:

$$\varepsilon^*(f) = \varepsilon_\infty + \frac{\Delta\varepsilon}{1 + j\,f/f_{\text{relax}}}$$

The transfer matrix for an E–L pair is computed, then raised to the N-pairs power via Chebyshev polynomial identity — **O(1) per frequency regardless of layer count**:

$$M_{\text{pair}}^P = U_{P-1}(s) \cdot M_{\text{pair}} - U_{P-2}(s) \cdot I$$

where $s = \frac{1}{2}\text{tr}(M_{\text{pair}})$ and $U_n$ are second-kind Chebyshev polynomials.

For a 120 mm thick sample with ~5 mm annual ring spacing, this represents **48 alternating layers**. The full derivation is in [`simulations/tmm_derivation.md`](simulations/tmm_derivation.md).

**Parameter reduction** — the key to well-posedness:

| Parameter | Strategy | Rationale |
|-----------|---------|-----------|
| $d_E$, $d_L$ | Fixed from annual ring measurement | Independent geometric constraint |
| $N_\text{pairs}$ | $N = \text{round}(D_\text{total} / (d_E + d_L))$ | Determined by total thickness |
| $\Delta\varepsilon_E$, $\Delta\varepsilon_L$ | Fixed (literature values) | Bounded by moisture content |
| $f_\text{relax}$ | **Shared** across E and L | Free-water relaxation dominates at 18–44 GHz |
| $\varepsilon_{\infty,E}$, $\varepsilon_{\infty,L}$ | **Optimized (2 params)** | Core unknowns — dry-wood matrix permittivity |

This yields a **3-parameter optimization**: $[\varepsilon_{\infty,E},\ \varepsilon_{\infty,L},\ f_{\text{relax}}]$.

### Two Polarization Modes

Under normal incidence, the orthotropic permittivity tensor $\bar{\bar{\varepsilon}}_r$ decouples into two independent scalar problems:

| Mode | E-field direction | Probes | Measurement |
|------|------------------|--------|------------|
| ∥ Parallel (∥) | Along fiber (L-axis) | $\varepsilon_L^*$ | S₂₁ with E ⊥ grain |
| ⟂ Perpendicular (⊥) | Across fiber (T-axis) | $\varepsilon_T^*$ | S₂₁ with E ∥ grain |

Modes are **measured separately and inverted independently**.

## Cost Function

### Pipeline 1 (power-law)

$$J = 0.5 \cdot \text{Huber}(|\text{S}_{21}|) + 0.25 \cdot \text{Huber}(\angle\text{S}_{21}) + 0.3 \cdot \tau_\text{err}^2 + 0.05 \cdot d_\text{reg}^2$$

- **Huber loss:** L2 for small residuals, L1 for large ones — robust to outlier frequency points
- **Group-delay constraint:** IFFT first-arrival peak couples optical path to thickness, breaking the d–tanδ degeneracy
- **Thickness regularization:** soft constraint within user-specified bounds

### Pipeline 2 (TMM + Debye)

Phase-enriched cost with Laplace likelihood (robust to outliers):

```matlab
% Magnitude: 0.5 * Huber loss (dB)
% Phase:     0.5 * Huber loss (radians)
```

MCMC posterior sampling uses `cost_tmm3_phase.m` (log-likelihood) + `log_prior_tmm3.m` (uniform prior with bounds).

## Optimization & Inference

### Deterministic Optimization

| Algorithm | Method | Typical pop |
|-----------|--------|-------------|
| Differential Evolution (DE) | Mutation F=0.5, Crossover CR=0.7 | 60 |
| CMA-ES | Covariance matrix adaptation | 10–100 |
| Hybrid DE-CMA-ES | Exchange elites every 60 generations | 100+100 |
| Genetic Algorithm (GA) | SBX crossover (η=21), mutation 0.15 | 30 |
| Particle Swarm (PSO) | ω=0.7, c₁=c₂=1.5 | 35 |

All use Latin Hypercube Sampling with Opposition-Based Learning (LHS-OBL) initialization.

### Bayesian MCMC

Adaptive Metropolis-Hastings (Haario et al., 2001):
- **4 parallel chains** — independent, randomly perturbed starts
- **5,000 burn-in** iterations with proposal covariance adaptation every 500 steps
- **20,000 post-burn-in samples** per chain (80,000 total posterior draws)
- **Split-$\hat{R}$** convergence diagnostic — all chains must achieve $\hat{R} < 1.1$
- **95% HDI** (Highest Density Interval) reported for all parameters

Posterior diagnostics (5 figures from `plot_posterior.m`):
1. Trace plots — 4 chains overlaid per parameter
2. Pairwise scatter (corner plot) — posterior correlations
3. Marginal histograms — per-parameter distributions
4. 2D density contours — bivariate posterior shape
5. Convergence diagnostics — $\hat{R}$ and acceptance rates

## Simulation Validation Suite

Six experiments validate the inversion pipeline before running on real data:

| # | Script | Purpose |
|---|--------|---------|
| Exp1 | `exp1_identifiability.m` | Cost landscape — confirm (tanδ, n) degeneracy |
| Exp2 | `exp2_recovery_test.m` | Monte Carlo recovery under 0.1–0.5 dB noise |
| Exp3 | `exp3_model_mismatch.m` | Fit Debye data with power-law model → systematic residuals |
| Exp4 | `exp4_debye_model.m` | Debye forward model validation |
| Exp5 | `exp5_fixed_d.m` | Well-posedness when thickness is fixed |
| Exp6 | `exp6_tmm_layered.m` | TMM vs uniform FP contrast sweep; 4-param reduced TMM recovery |

## Usage

### TMM + Debye pipeline (recommended)

```matlab
% 1. Run the main inversion (DE + MCMC)
run_tmm_inversion

% 2. Regenerate figures from saved results
%    (load results/TMM_Inversion_Results.mat or MCMC_Inversion_Results.mat)
group_meeting_figures
```

Outputs to `simulations/results/`:
- `TMM_Inversion_Fits.png` — measured vs fitted |S₂₁| for all 4 samples
- `TMM_Texture_Comparison.png` — vertical/parallel overlay per wood size
- `MCMC_Inversion_Results.mat` — full posterior chains + diagnostics
- 5 posterior diagnostic figures (if `plot_posterior` runs)

### Power-law pipeline (legacy)

```matlab
% Set WORK_MODE = 1 in run_wood_inversion.m, then:
run_wood_inversion     % batch inversion with 5 algorithms
plot_wood_comparison   % 4 comparison figures
n_scan_diagnosis       % diagnose frequency exponent n
```

### Simulation experiments

```matlab
cd simulations
exp1_identifiability    % → results/Exp1_*.png
exp2_recovery_test       % → results/Exp2_*.mat + .png
exp3_model_mismatch      % → results/Exp3_*.png
exp4_debye_model         % → results/Exp4_*.png
exp5_fixed_d             % → results/Exp5_*.png
exp6_tmm_layered         % → results/Exp6_*.png
synthesis                % print summary of all experiments
```

## Data Format

Keysight VNA CSV files. `load_measurement.m` auto-detects column names:

| Column | Expected names |
|--------|---------------|
| Frequency | `Freq_Hz_`, `Freq(Hz)`, `Freq`, `Frequency` |
| S₂₁ magnitude | `S21(DB)`, `S21_DB`, `S21 dB` |
| S₂₁ phase | `S21(DEG)`, `S21_DEG`, `Phase` |

Air reference files match `air*.csv`. Wood sample files are all other `.csv` files.

## Dependencies

- **MATLAB R2018b+**
- **Statistics and Machine Learning Toolbox** — `lhsdesign` (Latin hypercube sampling)
- **Signal Processing Toolbox** — `sgolayfilt` (Savitzky-Golay smoothing; code degrades gracefully if missing)
- **Parallel Computing Toolbox** — `parfor` (used in `run_tmm_inversion.m` for 4-sample parallel loop)
- **Optimization Toolbox** — `fmincon` (only needed by `n_scan_diagnosis.m`)

No third-party packages required.

## Key Research Findings

1. **4-parameter power-law model is underdetermined.** (tanδ, n) and (ε_r, d) are fundamentally coupled; |S₂₁| magnitude alone cannot separate them.

2. **Debye model with fixed thickness is well-posed.** Reducing from 4 free parameters to 3 (fixing d from caliper measurement, sharing f_relax across E/L layers) breaks the degeneracy.

3. **Relaxation frequency stable around 1–4 GHz**, consistent with bound-water Debye relaxation in wood (free-water relaxation at ~18 GHz is at the lower edge of the measurement band).

4. **Clear dielectric anisotropy:** ε_parallel systematically higher than ε_perpendicular, consistent with the known dielectric anisotropy of cellulose microfibrils.

5. **Layered structure matters only at high contrast.** exp6 shows that when Δε' between earlywood and latewood exceeds ~2.0, the TMM residual vs uniform FP exceeds the 0.1 dB detection threshold. For typical temperate wood with Δε' ≈ 1.5–2.0, the layered effect is subtle but real.

6. **Bayesian MCMC provides uncertainty quantification.** Posterior HDI widths give rigorous confidence bounds on inverted parameters, enabling statistical comparison between grain orientations.

## Hardcoded Paths

`run_wood_inversion.m` line 16–17 hardcodes `dataDir = 'D:\woods\lzmwoods'`. Several other scripts derive paths from the repo root using `pwd`. If you relocate the repo, update these paths or ensure you `cd` to the repo root before running scripts.

## License

MIT
