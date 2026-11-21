# Monte Carlo robustness analysis for a potential missing confounder

Reproducibility code for **A Monte Carlo procedure for quantifying the robustness of causal inference for a potential missing confounder**.

This repository contains the complete R implementation for the two illustrative examples, Figures 1 and 2, the 100-run statistics for Appendix A.5 Table A3, and the kindergarten-retention effect-size threshold sensitivity analysis for Appendix A.6 Table A4.

## Run the analysis

Only base R is required. The 2026-11-21 audit results were run with R 4.5.0 (2025-04-11), aarch64-apple-darwin20, on macOS 27.0.1. Analysis uses only base/recommended R packages: stats, graphics, grDevices, utils, datasets, methods and base (4.5.0); compiler/tools 4.5.0 appear in validation sessionInfo(). No contributed R packages are needed.

Download the repository or clone it, open a terminal in the repository directory, and run:

```sh
Rscript robustness.R
Rscript minimum_effect.R
Rscript retention_audit.R
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
| Education | 0.2984613894645406 | 0.22352057677494624 | 0.17686271729239786 | 10,562 | +1 |

## Reference results for Appendix Table A3

| Scenario | Mean admissible pairs | Mean probability | SD across 100 runs | MCSE of 100-run mean | First 50 mean | Last 50 mean |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Retention U(-1,1) | 594284.2 | 0.678885694346 | 0.00872590642718 | 0.000872590642718 | 0.679302751842 | 0.678468636851 |
| Retention U(-0.5,0.5) | 1000000 | 0.83659422 | 0.00652365187216 | 0.000652365187216 | 0.83601426 | 0.83717418 |
| Education U(-1,1) | 599629.39 | 0.825283992413 | 0.00781046584493 | 0.000781046584493 | 0.825522876062 | 0.825045108765 |
| Education U(-0.5,0.5) | 1000000 | 0.99697109 | 0.000397324285412 | 3.97324285412e-05 | 0.99692492 | 0.99701726 |

The final estimate is the arithmetic mean of 100 per-run conditional proportions. SD measures variability across complete runs. MCSE of the 100-run mean is SD/sqrt(100), not SD. Cross-combined pairs share marginal draws and are not independent Bernoulli observations; admissible-pair counts are not effective independent sample sizes. First/last-half comparisons are descriptive stability diagnostics, not a proof of convergence.

Figure labels use three decimal places: 0.679, 0.837, 0.825 and 0.997. Education changed from 0.820084 / 0.995132 to 0.825283992413386 / 0.99697109 after correcting the multi-control cross term. The restricted-prior result is below one; rounding is never evidence of exactly 100%. The earlier restricted-prior retention value 0.847420 was a single run; its repeated-run mean remains 0.83659422.

## Education calculation correction

`inputs.R` computes from Frank (2000), Table 1: `rxy=.427`, `rXZ=c(-.437,-.289)`, `rYZ=c(-.317,-.350)` and `RZZ=matrix(c(1,.263,.263,1),2)`. The cross term is `crossprod(rXZ, solve(RZZ,rYZ))`, **not** the product of the two multiple correlations. The squared multiple correlations are the corresponding quadratic forms. Thus:

```text
cross = 0.18838965505016483
rxz2 = 0.22352057677494624
ryz2 = 0.17686271729239786
partial_r = (rxy-cross)/sqrt((1-rxz2)*(1-ryz2)) = 0.2984613894645406
T(0,0) = partial_r*sqrt(10562/(1-partial_r^2)) = 32.138154101648
```

The text may display 0.30; calculations and figures use unrounded matrix results. Independent checks use a positive-definite full correlation matrix, its precision matrix, regression residuals and direct OLS. `inputs.R` also centralizes unchanged Retention inputs/scales for A.5 and A.6.


## Appendix A.6 effect-size threshold replacement and Table A4

Run `Rscript minimum_effect.R` (or `Rscript minimum_effect.R output_directory`) to generate `minimum_effect_runs.csv` and `minimum_effect_summary.csv`. This retention-only analysis replaces the original PMC critical value `d_alpha = 1.96` with the dimensionless, scenario-specific threshold `d_delta = Delta_min / SE_f`. The original NHST calculation remains the comparison baseline.

The supplied source for the reading-score standard deviation is Hong & Raudenbush (2005), Table 2 (p. 210), spring treatment-year reading score: `SD(Y) = 13.53`. Motivated by [Kraft (2020)](https://doi.org/10.3102/0013189X20912798), the code computes `delta_min <- 0.20 * 13.53`, giving an absolute threshold of 2.706 reading-score points without pre-rounding. The working scales remain `s_x = 0.240548` and `s_y = 13.53`. The published 13.53 is for the broader reading sample (N=11,797), not an identified at-risk or mixed-model residual SD; here it is an illustrative working scale. The threshold is divided by the full-model standard error separately for each admissible confounding pair; it is not inserted directly into the standardized PMC formula.

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

The six-step design, 1,000 draws per marginal, 100 repetitions, 1.96, 0.8, two uniform priors, and original RNG ordering are retained. Education inputs are now computed without pre-rounding. Mersenne-Twister / Inversion / Rejection uses seed 123: reconstruction retains its initial discarded `mcdat()` call; the four PMC scenarios reset seed 123 once and proceed in Retention-full, Retention-restricted, Education-full, Education-restricted order. A.6 resets once and continues between its priors.

`robustness_runs.csv` saves every repetition and finite-value diagnostics. All admissible configurations are checked, without silently dropping nonfinite values. Two complete runs produced byte-identical robustness CSVs and PNGs; A.6 CSVs and standalone diagnostic CSVs were also byte identical. PDF metadata timestamps can differ. `reproducibility.txt` records sessionInfo and validation commands/results.

## Retention source and covariance audit

See [SOURCE_AUDIT.md](SOURCE_AUDIT.md), [AUDIT_REPORT.md](AUDIT_REPORT.md), `retention_source_cells.csv` (all 30 source cells including the empty promoted stratum 15), `retention_cell_source_audit.csv` and `retention_diagnostics.csv`.

All 29 nonempty Table 8 logit count/mean/SD inputs match the old generator, and all Table 9 reading means/SDs match. The old generator combines Table 8 counts (471 retained + 7,168 promoted) with Table 9 outcome summaries (469 + 7,139). The 31-person mismatch affects nine cells. The reference is preserved as a synthetic full-count scenario, requiring missing readings to share the observed cell moments. The outcome-matched scenario uses N=7,608, but still assumes the Table 8 logit moments apply to the reading-observed subset. Neither reconstructs raw data.

`retention_reconstruction.R` exposes `mcdat(rho_within=0, sample="reference", return_data=FALSE)`. `sample="outcome_matched"` changes counts; a configurable common rho controls the unidentified within-cell covariance through correlated normal residuals. Published marginal means/SDs are distribution parameters; finite synthetic cells need not match them exactly. Zero rho preserves original draws and RNG state bit for bit. The 1–29 `strata` codes index cells; they are not the original 15 shared strata and must not be used as controls that encode treatment.

`Rscript retention_audit.R` uses 1,000 synthetic datasets per count/rho scenario and seed 123, resetting per scenario for common latent normals within each count design. Its fixed rho grid is -0.50, -0.25, 0, 0.25, 0.50. This is separate from the 100-repetition PMC calculation. The audit has no initial discarded simulation; hence its reference means differ slightly from the legacy 1,000-simulation stream. CSVs distinguish mean absolute discrepancy from absolute discrepancy of mean T, and include across-dataset SD, MCSE=SD/sqrt(1000), and analytic mixture moments.

Formal inputs stay `partial_r=-.13`, `rxz2=.2704`, `ryz2=.2025`, `nu=7635`: T(0,0)=-11.456422589690, absolute discrepancy 1.813577410310 from published -13.27. Reconstruction diagnostics use augmented-model nu=N-4 (7635 or 7604); the reduced OLS y~x+z instead has N-3 df. Neither df reproduces the published mixed-model estimator. The original Equation 4 adjusts 14 shared-stratum indicators, continuous propensity logit, and assessment duration, with correlated school random intercept and retention slope. Missing joint moments, school membership and duration are not available from these marginal tables. A rho scenario near -13.27 does not identify true covariance or validate this approximation; no rho or partial r was fitted to that t.

The original article was audited through its [complete text mirror](https://www.academia.edu/29274731/Effects_of_Kindergarten_Retention_Policy_on_Childrens_Cognitive_Growth_in_Reading_and_Mathematics); the original typeset PDF could not be visually checked. This access limit and unresolved identification assumptions are retained in the source audit. No original individual-level data are distributed or claimed recovered.

## Full validation commands

From the repository root (repeat run writes to a separate copied source directory):

```sh
Rscript tests/calculation_audit.R
Rscript tests/retention_audit.R
Rscript tests/minimum_effect.R
Rscript robustness.R
mkdir -p ../robustness-run2
cp *.R ../robustness-run2/
Rscript ../robustness-run2/robustness.R
Rscript tests/robustness_results.R . ../robustness-run2
Rscript minimum_effect.R ../a6-run1
Rscript minimum_effect.R ../a6-run2
Rscript tests/minimum_effect_results.R ../a6-run1 ../a6-run2
Rscript retention_audit.R
```

Results tests independently replay all 400 PMC repetitions and all 200 A.6 repetitions, verify V, numerical finiteness, exceedance counts, means/SD/MCSE, and compare two complete runs. A.6 additionally uses FWL coefficient/SE and effect-scale exceedance events. Matrix/pointwise tests include boundaries and OLS checks. The per-run Cartesian pairs share marginal draws; no binomial MCSE is used.

Figure 1/2 design remains the existing black-and-white base-R design, PDF plus 600-dpi PNG. Render-time panel geometry is saved in `*.alignment-layout.json`; both panels pass the 1.5-point alignment test, rendered PDF text is at least 9 pt, and collision audits found zero failures/warnings. Generic Nature static checks request SVG/TIFF/Cairo devices; these format recommendations do not replace the explicitly requested existing PDF/PNG design. No journal-specific submission compliance is claimed.
