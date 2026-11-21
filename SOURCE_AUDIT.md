# Hong & Raudenbush (2005): source audit

Accessed 2026-11-21. Audited the original article, *Effects of Kindergarten Retention Policy on Children’s Cognitive Growth in Reading and Mathematics*, EEPA 27(3), 205–224, DOI [10.3102/01623737027003205](https://doi.org/10.3102/01623737027003205), using the [complete original article text on Academia](https://www.academia.edu/29274731/Effects_of_Kindergarten_Retention_Policy_on_Childrens_Cognitive_Growth_in_Reading_and_Mathematics). The publisher PDF route redirects to an abstract; the historical AERA PDF URL was inaccessible. Thus the transcription below is verified against accessible original-article text, not visually against the typeset PDF. The local repository `retention.pdf` is R-generated graphics, not the original article.

## Original-source facts

Both Tables 8 and 9 appear on printed page 216. Table 8 contains logit summaries; Table 9 contains reading summaries. Table 11, page 217, reports retention coefficient −9.01, SE 0.68, and t −13.27. Equation 4 includes treatment, continuous estimated propensity logit, 14 indicators for 15 strata, assessment duration, correlated school random intercept and treatment slope, and individual residual error. Table 11’s reported random-effect variances are 55.02 (intercept), 18.83 (slope), and 88.22 (individual); intercept–slope correlation is −0.27. Its note explicitly says the strata, logit and duration coefficients are omitted from the table. The individual propensity model uses 39 predictors and seven quadratic terms; 207 pretreatment covariates were assessed for balance. Eight retained children in stratum 15 have no promoted match. Table 2, printed page 210, gives spring treatment-year reading SD 13.53 for 11,797 children. [Original article text](https://www.academia.edu/29274731/Effects_of_Kindergarten_Retention_Policy_on_Childrens_Cognitive_Growth_in_Reading_and_Mathematics).

## Exact cell transcription and comparison

The CSV `retention_source_cells.csv` contains all 30 treatment-by-stratum cells, including the empty promoted stratum 15. Empty means/SDs are blank, not zero. All original-source means and SDs were entered separately from the R file, then compared programmatically.

Reference implementation: pre-correction `robustness.R`, commit `cd24f33`, `mcdat()`. The transcription is retained in `retention_source_cells.csv`.
- All 29 nonempty Table 8 counts, means and SDs match `mcdat()` exactly.
- All 29 nonempty Table 9 means and SDs match `mcdat()` exactly.
- `mcdat()` applies Table 8 counts to Table 9 reading moments. Nine counts disagree with the outcome table:

| Treatment | Stratum | Table 8 / mcdat N | Table 9 N | Deficit |
|---|---:|---:|---:|---:|
| retained | 9 | 46 | 45 | 1 |
| retained | 14 | 39 | 38 | 1 |
| promoted | 1 | 3,054 | 3,044 | 10 |
| promoted | 2 | 1,661 | 1,654 | 7 |
| promoted | 3 | 983 | 977 | 6 |
| promoted | 4 | 323 | 321 | 2 |
| promoted | 5 | 441 | 440 | 1 |
| promoted | 6 | 154 | 153 | 1 |
| promoted | 7 | 213 | 211 | 2 |

All remaining cells’ counts agree. Table 8 totals are 471 retained and 7,168 promoted, N=7,639. Table 9 totals are 469 retained and 7,139 promoted, N=7,608. The difference is 31, comprising two retained and 29 promoted. These totals and differences are independently computed from the transcribed cells.

## Implications and identification boundaries (audit analysis)

The original 7,639-count generator is a synthetic full-count scenario. It effectively extends observed reading-cell means/SDs to the 31 additional observations. That needs a missing-outcome assumption: conditional on treatment and stratum, the unobserved outcomes must share the supplied first two moments if those moments are to characterize the full cells. Conditional MCAR or an appropriate MAR model could motivate such a scenario, but the count discrepancy alone does not establish either. Missing outcomes could have different moments.

A 7,608-count scenario respects Table 9’s observed-outcome counts. However, Table 8 does not give the propensity-logit means/SDs specifically among those with observed reading outcomes. Reusing Table 8 logit moments after changing counts therefore remains an additional approximation. It cannot be described as recovering an exact complete-case correlation matrix.

Neither table supplies within-cell covariance of reading and propensity logit. Drawing z and y independently within each cell assumes that covariance is zero. This is a modeling choice, not a consequence of the displayed marginal summaries. Their Gaussian generation also does not recover the original child-level distribution. Global covariance can still arise from differences among cell means even with zero within-cell covariance. Varying this unidentified covariance changes reconstructed marginal/partial correlations while preserving the supplied within-cell marginals.

The two treatment groups must share the same 15 strata. The current generator assigns distinct `strata` labels 1–29 to treatment-by-stratum cells; those labels encode treatment group as well as stratum. They may index cells, but are not equivalent to the article’s 15 shared propensity strata. Using all 29 cell indicators as controls makes treatment collinear with the controls.

Equation 4’s mixed model coefficient/SE cannot be recovered by ordinary least squares regression on treatment plus one continuous z, or by `r_xy.z`. Adding correct shared stratum indicators would still omit duration, school membership, school random-effect structure and unknown covariance. These missing ingredients are not identified by the two tables. A conversion t/sqrt(t²+nu) is an ordinary-regression identity with the corresponding residual df; assigning nu=7,635 does not turn the article’s mixed-model t into the partial correlation of this synthetic three-variable regression. For OLS with an intercept, treatment and one continuous z, residual df would be N−3 (7,636 or 7,605). N−4 (7,635 or 7,604) corresponds instead to a model with an additional scalar predictor, such as the hypothetical confounder. Neither establishes the correct df for Equation 4.

The SD 13.53 describes the broader spring reading sample in Table 2. It is not established by that source as the reading SD among 7,639 at-risk children, the observed-outcome 7,608, a within-cell SD, or a mixed-model residual SD. Reusing it needs to be labeled as a chosen working scale. Likewise, the 8 unmatched children in stratum 15 do not establish overlap; their inclusion in a simulation cannot repair that support limitation.

## Source limitations

No raw data, reading-observed propensity moments, within-cell cross moments, school identifiers, assessment times or complete Equation 4 fitted design/covariance matrix were obtained. No claim of exact study replication is supported. A different 2006 JASA paper, *Evaluating Kindergarten Retention Policy*, is not interchangeable with this 2005 EEPA article and was not used for the numerical transcription.
