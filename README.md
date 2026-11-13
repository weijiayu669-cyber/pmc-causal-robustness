# Monte Carlo robustness analysis for a potential missing confounder

Reproducibility code for **A Monte Carlo procedure for quantifying the robustness of causal inference for a potential missing confounder**.

This repository contains the complete R implementation for the two illustrative examples, Figures 1 and 2, the 100-run statistics for Appendix A.5 Table A3, and the kindergarten-retention effect-size threshold sensitivity analysis for Appendix A.6 Table A4.

## Run the analysis

Only base R is required. The reference results were verified with R 4.5.0 on macOS.

Download the repository or clone it, open a terminal in the repository directory, and run:

```sh
Rscript robustness.R
Rscript minimum_effect.R
```

Alternatively, in R or RStudio:

```r
source("robustness.R")
source("minimum_effect.R")
```

Run `robustness.R` first to reproduce Table A3 and Figures 1 and 2, then `minimum_effect.R` to reproduce Table A4. The scripts are standalone and do not depend on each other's output files.

`robustness.R` writes the original summary CSV and both figures beside itself, replacing the included reference outputs. It prints correlation reconstruction and numerical stability diagnostics to the console. Its complete run took about one minute on the reference machine; runtime varies by hardware. `minimum_effect.R` writes the Table A4 CSV and the per-run A.6 CSV listed below beside itself; it does not write the original Table A3 CSV or figures.

## Files

| File | Contents |
| --- | --- |
| `robustness.R` | Standalone computation, 1,000 correlation simulations, four 100-run scenarios, and figure generation |
| `minimum_effect.R` | Standalone Appendix A.6 effect-size threshold analysis for kindergarten retention, 100 runs per prior |
| `robustness_summary.csv` | Appendix A.5 Table A3 estimates, SD, MCSE, admissible-pair counts, and first/last-half diagnostics |
| `minimum_effect_summary.csv` | Appendix A.6 Table A4: original NHST and effect-size threshold scenario proportions, SD, MCSE, admissible counts and percentage-point differences |
| `minimum_effect_runs.csv` | Appendix A.6: all 100 runs per prior, paired baseline/effect-size threshold results and numerical diagnostics |
| `README.md` | Running order, R version, random seeds, model scope, and correspondence between scripts and appendix tables |
| `retention.pdf` / `.png` | Kindergarten retention, full and restricted priors |
| `education.pdf` / `.png` | Educational attainment, full and restricted priors |

PDF figures use vector contours and text; PNG figures are 7,200 by 4,200 pixels (12 by 7 inches at 600 DPI). Both retain the original monochrome contours, shading, and legend.

## Implemented calculation

1. Obtain the observed partial correlation, squared covariate correlations, original effect direction, and residual degrees of freedom. The two supplied examples use the original working inputs listed below.
2. Set the critical value to 1.96 and the pointwise power threshold to 0.8.
3. Specify independent uniform nominal marginals on either [-1, 1] or [-0.5, 0.5].
4. Draw 1,000 values for each marginal and cross-combine them. Retain pairs satisfying both positive residual-variance conditions and a positive determinant.
5. Calculate directional pointwise power for each admissible pair and the proportion whose power exceeds 0.8.
6. Interpret the result conditional on admissibility and the specified prior. Repeat the complete simulation 100 times per scenario to assess numerical precision, and generate the two contour figures.

For a pair `(x, y)` of hypothetical confounder correlations, the shared calculation is:

```text
q = partial_r * sqrt((1 - rxz2) * (1 - ryz2)) - x*y
vx = 1 - rxz2 - x*x
vy = 1 - ryz2 - y*y
D = vx*vy - q*q
admissible: vx > 0, vy > 0, D > 0
p_mc = Phi(direction * q / sqrt(D / nu) - 1.96)
```

`Phi` denotes the standard normal cumulative distribution function. The scenario probability is `Pr(p_mc > 0.8 | admissible)`; it is distinct from the mean pointwise power.

| Example | Partial correlation | rxz2 | ryz2 | nu | Direction |
| --- | ---: | ---: | ---: | ---: | ---: |
| Retention | -0.13 | 0.2704 | 0.2025 | 7,635 | -1 |
| Education | 0.29 | 0.224 | 0.177 | 10,562 | +1 |

## Reference results for Appendix Table A3

| Scenario | Mean admissible pairs | Mean probability | SD across 100 runs | MCSE of 100-run mean | First 50 mean | Last 50 mean |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Retention U(-1,1) | 594,284.2 | 0.678886 | 0.008725906 | 0.000872591 | 0.679303 | 0.678469 |
| Retention U(-0.5,0.5) | 1,000,000 | 0.836594 | 0.006523652 | 0.000652365 | 0.836014 | 0.837174 |
| Education U(-1,1) | 601,030.7 | 0.820084 | 0.007962974 | 0.000796297 | 0.820310 | 0.819859 |
| Education U(-0.5,0.5) | 1,000,000 | 0.995132 | 0.000579725 | 0.000057973 | 0.995075 | 0.995189 |

The final estimate is the arithmetic mean of 100 per-run conditional proportions. SD measures variability across complete runs. MCSE of the 100-run mean is SD/sqrt(100), not SD. Cross-combined pairs share marginal draws and are not independent Bernoulli observations; admissible-pair counts are not effective independent sample sizes. First/last-half comparisons are descriptive stability diagnostics, not a proof of convergence.

Figure labels use these 100-run means rounded to two decimal places: 0.68, 0.84, 0.82, and 1.00. The last value is rounded from 0.995132 and does not mean exactly 100%. The earlier restricted-prior retention value 0.847420 was a single-run estimate; the final repeated-run mean is 0.83659422.

## Appendix A.6 effect-size threshold replacement and Table A4

Run `Rscript minimum_effect.R` (or `Rscript minimum_effect.R output_directory`) to generate `minimum_effect_runs.csv` and `minimum_effect_summary.csv`. This retention-only analysis replaces the original PMC critical value `d_alpha = 1.96` with the dimensionless, scenario-specific threshold `d_delta = Delta_min / SE_f`. The original NHST calculation remains the comparison baseline.

The supplied source for the reading-score standard deviation is Hong & Raudenbush (2005), Table 2A, spring treatment-year reading score: `SD(Y) = 13.53`. Motivated by [Kraft (2020)](https://doi.org/10.3102/0013189X20912798), the code computes `delta_min <- 0.20 * 13.53`, giving an absolute threshold of 2.706 reading-score points without pre-rounding. The working scales are `s_x = 0.240548` and `s_y = 13.53`. The threshold is divided by the full-model standard error separately for each admissible confounding pair; it is not inserted directly into the standardized PMC formula.

With `a = r_xc`, `b = r_yc`, and the existing assumption that the omitted confounder is orthogonal to observed covariates:

```text
q = -0.13 * sqrt((1 - 0.52^2) * (1 - (-0.45)^2))
vx = 1 - 0.52^2 - a^2
vy = 1 - (-0.45)^2 - b^2
A = q - a*b
D = vx*vy - A^2
admissible V: vx > 0, vy > 0, D > 0
tt = A / sqrt(D / 7635)
SE_f = (13.53 / 0.240548) * sqrt(D / 7635) / vx
Delta_min = 0.20 * 13.53
d_delta = Delta_min / SE_f
p0 = Phi(-tt - 1.96)
p_delta = Phi(-tt - d_delta)
Q0 = Pr(p0 > 0.8 | V)
Q_delta = Pr(p_delta > 0.8 | V)
```

All divisions and square roots for `tt`, `SE_f` and the probabilities are evaluated only on admissible pairs. Equation (A10) is consistent with the original statistic: `beta_f = (s_y/s_x)*A/vx` implies `beta_f/SE_f = tt = sqrt(7635)*A/sqrt(D)`. The production evaluator does not need to construct `beta_f`. The reduced-model effect direction remains negative.

The replacement is exactly `1.96 -> Delta_min/SE_f`; no extra `-1.96` or additional significance condition is applied. Both `SE_f` and `d_delta` vary with `(a,b)`. Unlike a constant cutoff of 2.706, the new probability has a coherent effect-scale interpretation under the specified normal approximation: for `B ~ N(beta_f, SE_f^2)`, it equals `Pr(B < -Delta_min)`, holding the scenario-specific standard error fixed. This is not a posterior probability about the true causal effect. No universal ordering between `p_delta` and `p0` is imposed, since `Delta_min/SE_f` may be above or below 1.96.

The Monte Carlo design is unchanged: independent nominal uniform marginals on `[-1,1]` or `[-0.5,0.5]`, 1,000 draws per correlation, 1,000,000 Cartesian candidate pairs per run, and 100 runs per prior. Only admissible pairs enter the proportions. Mersenne-Twister / Inversion / Rejection and seed 123 are set once before the full prior; the restricted prior continues the same stream without resetting, matching Appendix A.5.

| Prior | Original Q0 | Effect-size Q_delta | Across-run SD of Q_delta | MCSE | Mean admissible pairs | Q0 - Q_delta (percentage points) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| U(-1,1) | 0.678885694346454 | 0.647879036365835 | 0.00885235356829177 | 0.000885235356829177 | 594284.2 | 3.1006657980619 |
| U(-0.5,0.5) | 0.83659422 | 0.79485314 | 0.00694559479315538 | 0.000694559479315538 | 1e+06 | 4.174108 |

`minimum_effect_runs.csv` contains 200 rows (100 per prior): both proportions, exceedance counts, admissible counts, paired percentage-point differences, fixed scale inputs and finite-value diagnostics. Per-run minimum/maximum values of `SE_f` and `d_delta` document their variation. `minimum_effect_summary.csv` contains the two Table A4 rows plus baseline SD/MCSE, paired-difference uncertainty and first/last-half diagnostics. Both CSVs retain approximately 15 significant digits. MCSE is the across-run SD divided by `sqrt(100)`; shared marginal draws mean the Cartesian pairs are not independent Bernoulli observations. Mean pointwise probabilities are recorded separately from the scenario proportions.

Validation commands, from the repository root:

```sh
Rscript tests/minimum_effect.R
Rscript minimum_effect.R /tmp/pmc-scaled-run1
Rscript minimum_effect.R /tmp/pmc-scaled-run2
Rscript tests/minimum_effect_results.R /tmp/pmc-scaled-run1 /tmp/pmc-scaled-run2
```

The pointwise tests compare `SE_f` and the unchanged T-statistic with direct `lm()` fits having exactly specified sample covariance matrices. They recover NHST when `Delta_min = 1.96*SE_f`, verify probability 0.5 at `beta_f = -Delta_min` (detecting an extra significance cutoff), and verify invariance when the reading-score scale and absolute threshold are rescaled together. Full correlation-matrix positive definiteness independently checks admissibility, with separate boundary tests. Nonfinite values on admissible pairs cause failure rather than silent removal.

The results test verifies byte-identical CSVs from two complete runs, summary arithmetic, Appendix A.5 statistics and all 200 baseline runs against the unchanged main evaluator. It independently reconstructs the coefficient and SE using FWL residual cross-products and verifies every new proportion through the equivalent event `beta_f + qnorm(0.8)*SE_f < -2.706`. First/last-half comparisons are descriptive stability diagnostics, not a convergence proof.

The original Q0 means reproduce Appendix A.5 at full stored precision and round to 0.678886 and 0.836594. The new Table A4 values are recomputed with `s_y = 13.53` and `Delta_min = 2.706`; they supersede earlier constant-cutoff and differently scaled results. This analysis retains the original working OLS model and does not refit the source study's multilevel estimator.

## Correlation simulation for the retention example

The script generates 1,000 synthetic datasets from the embedded stratum-specific summaries. The verified correlation means are:

| Correlation | Mean | First 500 mean | Last 500 mean |
| --- | ---: | ---: | ---: |
| r_xz | 0.523155 | 0.523203 | 0.523108 |
| r_yz | -0.450598 | -0.450678 | -0.450519 |
| r_xy | -0.332705 | -0.332838 | -0.332573 |
| r_xy given z | -0.127458 | -0.127560 | -0.127357 |

Within each treatment-by-stratum cell, propensity-score logits and reading scores are generated independently from the specified normal distributions. These are assumption-dependent working correlations, not a uniquely identified original-study correlation matrix. The rounded inputs 0.52, -0.45, and -0.13 are retained in the PMC calculation; the covariate correlations enter it as squares.

## Reproducibility details

The script preserves the mathematical formulas, rounded working inputs, degrees of freedom, and random-number order of the supplied implementations. It uses Mersenne-Twister / Inversion / Rejection and two seed-123 streams:

- The correlation reconstruction and original single-run comparisons follow the original figure script, including its initial `mcdat()` call before the 1,000 recorded simulations.
- The four 100-run scenarios restart from `set.seed(123)`, matching the original standalone diagnostics script.

The merged implementation was checked against all 1,000 recorded correlation simulations and all 400 scenario repetitions. A separate standalone run reproduced the summary CSV and PNG files exactly on the reference machine. Rendering may differ across operating systems or graphics devices.

All numerical inputs needed to run the examples are embedded in the script; no individual-level study data are distributed. The stratum-specific summary inputs are retained as supplied. Their exact table numbering and the reading-outcome versus propensity-score sample counts require confirmation against the final source publication; this repository does not change the original N=7,639 or residual degrees of freedom. Numerical repeatability does not establish equivalence to the source study's multilevel estimator.
