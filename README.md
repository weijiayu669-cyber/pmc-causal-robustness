# Monte Carlo robustness analysis for a potential missing confounder

Reproducibility code for **A Monte Carlo procedure for quantifying the robustness of causal inference for a potential missing confounder**.

This repository contains the complete R implementation for the two illustrative examples, Figures 1 and 2, the 100-run statistics for Appendix A.5 Table A3, and the kindergarten-retention minimum-effect-size sensitivity analysis for Appendix A.6 Table A4.

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

`robustness.R` writes the original summary CSV and both figures beside itself, replacing the included reference outputs. It prints correlation reconstruction and numerical stability diagnostics to the console. Its complete run took about one minute on the reference machine; runtime varies by hardware. `minimum_effect.R` writes the Table A4 CSV and the additional A.6 diagnostics listed below beside itself; it does not write the original Table A3 CSV or figures.

## Files

| File | Contents |
| --- | --- |
| `robustness.R` | Standalone computation, 1,000 correlation simulations, four 100-run scenarios, and figure generation |
| `minimum_effect.R` | Standalone Appendix A.6 minimum-effect analysis for kindergarten retention, 100 runs per prior |
| `robustness_summary.csv` | Appendix A.5 Table A3 estimates, SD, MCSE, admissible-pair counts, and first/last-half diagnostics |
| `minimum_effect_summary.csv` | Appendix A.6 Table A4: original NHST and minimum-effect directional-power scenario proportions, SD, MCSE, and admissible counts |
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

## Appendix A.6 minimum-effect-size analysis and Table A4

Run `Rscript minimum_effect.R` to reproduce the included `minimum_effect_summary.csv`. This extension is restricted to kindergarten retention and preserves the original working correlations, residual degrees of freedom, admissible region, cutoff 1.96, and two uniform nominal priors. Each run draws 1,000 values for each correlation and forms their Cartesian product; the calculation is repeated 100 times per prior. It uses Mersenne-Twister / Inversion / Rejection, resets seed 123 once before the full prior, then runs the restricted prior without another reset, matching the retention portion of the original diagnostics stream.

The illustrative minimum-effect threshold is `Delta = 0.20 * (9.01 / 0.67) = 2.689552...` reading-score points. The 0.67-SD description is approximate. Motivated by [Kraft (2020)](https://doi.org/10.3102/0013189X20912798), 0.20 SD is an analyst-selected sensitivity threshold; it is not a policy consensus or a universal minimum-effect standard.

For clarity, the following notation uses `q` for the residual covariance before the omitted-confounder adjustment, and `A` for the adjusted numerator (the latter is called `q` in the original calculation above):

```text
q = -0.13 * sqrt((1 - 0.52^2) * (1 - 0.45^2))
vx = 1 - 0.52^2 - a^2
vy = 1 - 0.45^2 - b^2
A = q - a*b
D = vx*vy - A^2
K = s_y / s_x
beta_c = K * A / vx
SE_c = K * sqrt(D / 7635) / vx
p_mc = Phi(-beta_c / SE_c - 1.96)
p_delta = Phi((-Delta - beta_c) / SE_c - 1.96)
```

These coefficient and conventional OLS standard-error formulas follow from the Schur complement and Frisch-Waugh-Lovell regression under the existing assumption that the omitted confounder is orthogonal to the observed covariates. Their ratio is exactly `A / sqrt(D / 7635)`, so the original pointwise power is unchanged. The minimum-effect test is `H0: beta_c >= -Delta` versus `H1: beta_c < -Delta`; retaining the 1.96 cutoff gives a directional normal-test size of approximately 0.025. The code holds the scenario-specific standard error fixed in this normal power approximation.

The score scale is fixed from the supplied code's group sizes, reading means and within-group SDs. With `N = 7639`, the pooled-summary identity `s_y^2 = sum((n_g - 1)*s_g^2 + n_g*(mean_g - mean_y)^2)/(N - 1)` gives `s_y = 13.374335...`; binary retention gives `s_x = sqrt(N/(N - 1)*(471/N)*(7168/N)) = 0.240547789...`. Thus `K = 55.599493...`. This is an explicitly defined working scale, not a recovered exact source-study SD. The original study's multilevel coefficient and standard error cannot be identified from the three working correlations and these summaries. This extension does not refit that multilevel model or force its working coefficient to equal -9.01.

The primary summary is `Pr(p_delta > 0.8 | admissible)`, averaged over the 100 per-run proportions. It is distinct from mean pointwise power and from filtering scenarios by a plug-in effect estimate.

| Prior | Criterion | Mean scenario proportion | SD across runs | MCSE | Mean admissible count |
| --- | --- | ---: | ---: | ---: | ---: |
| U(-1,1) | Original NHST directional power | 0.678886 | 0.008726 | 0.000873 | 594,284.2 |
| U(-1,1) | Minimum-effect directional power | 0.616517 | 0.008681 | 0.000868 | 594,284.2 |
| U(-0.5,0.5) | Original NHST directional power | 0.836594 | 0.006524 | 0.000652 | 1,000,000 |
| U(-0.5,0.5) | Minimum-effect directional power | 0.744388 | 0.007044 | 0.000704 | 1,000,000 |

`minimum_effect_summary.csv` contains these four Table A4 rows, with full numerical precision and additional diagnostic columns. `mean` is the mean scenario proportion; `mean_pointwise_power` is a separate quantity. MCSE equals the SD of 100 complete runs divided by 10. Within a run, the same configurations are reused across criteria; shared marginal draws mean the million candidate pairs are not independent observations.

The script also regenerates these supporting files locally; they need not be downloaded separately:

- `A6_runs_100.csv`: 100 runs per prior, with five criteria recorded on each run's common configurations (1,000 long-format rows, not 1,000 independent runs).
- `A6_summary.csv`: all five criteria, including the original power, the primary minimum-effect power, a joint sampling-event power, a plug-in effect filter, and a working-SD threshold check.
- `A5_baseline_reproduction.csv`: retention baseline comparison with Table A3, including precision and first/last-half checks.
- `A6_comparison.csv`: new-minus-original paired differences and their run-based SD and MCSE.
- `scale_audit.csv`, `A6_sampled_axes.rds`, and `sessionInfo.txt`: scale inputs, all sampled correlation axes, and the R environment.

For the secondary joint sampling event, `B ~ N(beta_c, SE_c^2)` must satisfy both `B < -1.96*SE_c` and `B < -Delta`; its power is `Phi(-beta_c/SE_c - max(1.96, Delta/SE_c))`. The resulting scenario proportions are 0.647374 and 0.794357. The deterministic filter `p_mc > 0.8` and `-beta_c >= Delta` gives 0.658703 and 0.812550. Neither is the primary minimum-effect test. Using exactly 0.20 of the working SD instead gives `Delta = 2.674867...` points and primary proportions of 0.616914 and 0.745022; this is explicitly a different SD basis.

The A.6 baseline was checked against the original evaluator for all 200 retention runs, and the coefficient/standard-error derivation was checked against direct `lm()` fits with exactly specified covariance matrices. The zero-threshold identity and pointwise ordering `p_delta <= p_joint <= p_mc` were also verified. Numerical repeatability validates the implementation within its working model; it does not supply the individual-level design, school covariance, weights, or confounder associations needed for an exact reanalysis of the original multilevel estimator.

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
