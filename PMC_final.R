# PMC_final.R -- base R only; run with Rscript PMC_final.R or source().
# Outputs beside this script: one Table A3 CSV and Figures 1/2 as PDF + PNG.
# Mathematical inputs, degrees of freedom, 1.96 cutoff, 0.8 threshold,
# 1,000 marginal draws, 100 repetitions and both seed-123 streams are retained.
# Correlation reconstruction and original single-run checks use final_ver's RNG
# order, INCLUDING its initial mcdat() call before the 1,000 replications.
# Table A3 then resets seed=123, exactly as the separate diagnostics script.
# Figure design is unchanged; probability labels now use Table A3's 100-run means.
# The four plotted labels are rounded to two decimals, as in final_ver.

output_dir <- local({
  source_files <- lapply(sys.frames(), function(frame) frame$ofile)
  source_files <- Filter(function(path) is.character(path) && length(path) == 1L,
                         source_files)
  script_arg <- grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
  script_file <- if (length(source_files)) {
    tail(source_files, 1L)[[1L]]
  } else if (length(script_arg)) {
    sub('^--file=', '', script_arg[1L])
  } else {
    NULL
  }
  script_dir <- if (is.null(script_file)) getwd() else dirname(normalizePath(script_file))
  path <- script_dir
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
})

vrnorm <- Vectorize(rnorm,c("n","mean","sd"))

mcdat <- function() {
  dat <- mat.or.vec(7639,4)
  size <- c(9,12,14,12,24,23,47,48,46,48,47,49,45,39,8,3054,1661,983,323,441,154,213,143,85,49,35,16,8,3)
  dat[,1] <- rep(c(1,0),times=c(471,7168))
  dat[,2] <- rep(1:29,times=size)
  # Define means and standard deviations for the logit of propensity scores (z)
  meanz <- c(-5.42, -4.25, -3.42, -2.98, -2.57, -2.15, -1.75, -1.24, -0.83, -0.46, 
             -0.01, 0.53, 1.16, 1.99, 3.7, -5.38, -4.27, -3.5, -2.97, -2.55, 
             -2.14, -1.77, -1.27, -0.86, -0.47, -0.04, 0.56, 1.14, 1.8)
  sdz <- c(0.4, 0.24, 0.2, 0.1, 0.17, 0.07, 0.15, 0.14, 0.11, 0.12, 0.11, 0.22, 
           0.2, 0.33, 1.09, 0.39, 0.25, 0.2, 0.1, 0.16, 0.07, 0.15, 0.14, 0.1, 
           0.11, 0.15, 0.17, 0.17, 0.17)
  z <- vrnorm(size,meanz,sdz)
  z <- unlist(z)
  dat[,3] <- z
  # Define means and standard deviations for reading scores (y)
  meany <- c(47.41,47.25,42.66, 45.24,44.07,36.63,37.41,39.79,37.78,35,34.3,32.06,31.09,33.87,32.84,
             59.22,55.62,52.49,52.45,50.05,48.77,45.89,46,43.28,42.18,40.45,33.03,35.27,40.47)
  sdy <- c(11.73,15.50,15.92,17.32,14.43,13.88,10.92,15.37,12.27,9.77,9.05,9.31,11.96,9.15,
           7.71,11.64,11.30,11.95,12.53,12.15,12.67,12.40,13.37,11.52,13.34,13.08,14.33,8.25,8.60)
  y <- vrnorm(size,meany,sdy)
  y <- unlist(y)
  dat[,4] <- y
  dat <- data.frame(dat)
  colnames(dat) <- c("x","strata","z","y")
  rxz <- cor(dat$x,dat$z)
  ryz <- cor(dat$y,dat$z)
  rxy <- cor(dat$x,dat$y)
  rxyonz <- (rxy-rxz*ryz)/sqrt((1-rxz^2)*(1-ryz^2))
  result <- c(rxz,ryz,rxy,rxyonz)
  return(result)
}

# Explicitly retain the default RNG used in the supplied R scripts (R >= 3.6).
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(123)
cat("R version:", R.version.string, "\nRNG kind:",
    paste(RNGkind(), collapse=", "), "\nSeeds: 123 for each original script stream.\n")

# The first 15 group sizes total 471 retained children; the other 14 total
# 7,168 promoted children (N=7,639). Keep all original summary inputs.
# z and y are generated independently WITHIN each treatment-by-stratum cell.
# This imposes zero within-cell covariance; Table 8's marginal summaries do
# not identify that covariance or a unique study correlation matrix.
# In particular, r_yz is negative. The downstream model uses r_yz^2.
# These are synthetic working correlations, not recovered raw study records.
initial_mcdat <- mcdat()  # Preserve original RNG consumption; not one of the 1,000.
output <- t(replicate(1000L, mcdat()))
colnames(output) <- c("rxz", "ryz", "rxy", "rxy|z")
correlation_summary <- data.frame(
  correlation=colnames(output), mean=colMeans(output),
  sd_across_simulations=apply(output, 2L, sd),
  mcse_mean_1000=apply(output, 2L, sd)/sqrt(1000L),
  first_500_mean=colMeans(output[1:500, ]),
  last_500_mean=colMeans(output[501:1000, ]), row.names=NULL)
correlation_summary$half_difference <- abs(correlation_summary$first_500_mean -
                                           correlation_summary$last_500_mean)
checkpoints <- c(100L, 250L, 500L, 750L, 1000L)
correlation_running_means <- data.frame(
  simulations=checkpoints,
  t(vapply(checkpoints, function(n) colMeans(output[seq_len(n), ]), numeric(4L))),
  check.names=FALSE)
cat("\nHong & Raudenbush: 1,000 synthetic data sets\n")
print(correlation_summary, digits=10, row.names=FALSE)
cat("\nCumulative means (stability diagnostic, not proof of convergence):\n")
print(correlation_running_means, digits=10, row.names=FALSE)
cat("\nPMC working inputs retained: partial r=-0.13, r_xz=0.52, |r_yz|=0.45.\n",
    "Reconstruction r_yz is negative; squared inputs remain 0.2704 / 0.2025.\n",
    "Within-cell covariance is unidentified and is set to zero by mcdat().\n", sep="")

n_marginal <- 1000L
n_repeats <- 100L
threshold <- 0.8
scenarios <- list(
  list(name = "Retention U(-1,1)", bounds = c(-1,1), partial_r=-0.13,
       rxz2=0.52^2, ryz2=0.45^2, nu=7635, direction=-1),
  list(name = "Retention U(-0.5,0.5)", bounds = c(-0.5,0.5), partial_r=-0.13,
       rxz2=0.52^2, ryz2=0.45^2, nu=7635, direction=-1),
  list(name = "Education U(-1,1)", bounds = c(-1,1), partial_r=0.29,
       rxz2=0.224, ryz2=0.177, nu=10562, direction=+1),
  list(name = "Education U(-0.5,0.5)", bounds = c(-0.5,0.5), partial_r=0.29,
       rxz2=0.224, ryz2=0.177, nu=10562, direction=+1)
)

# Shared pointwise calculation, algebraically identical to hong()/edu().
# Strict admissibility is identical to diagnostics and the drawing surface.
# No binomial standard error is computed from Cartesian pairs: they share axes.
evaluate_pairs <- function(xs, ys, partial_r, rxz2, ryz2, nu, direction) {
  xx <- matrix(rep(xs, times=length(ys)), nrow=length(xs))
  yy <- matrix(rep(ys, each=length(xs)), nrow=length(xs))
  q <- partial_r * sqrt((1-rxz2) * (1-ryz2))
  vx <- 1 - rxz2 - xx*xx
  vy <- 1 - ryz2 - yy*yy
  adjusted <- q - xx*yy
  determinant <- vx*vy - adjusted*adjusted
  valid <- vx > 0 & vy > 0 & determinant > 0
  n_valid <- sum(valid)
  if (!n_valid) stop("No admissible candidate pairs")
  z <- direction * adjusted[valid] / sqrt(determinant[valid] / nu) - 1.96
  p <- pnorm(z)
  c(p_hat=mean(p > threshold), n_admissible=n_valid,
    acceptance_rate=n_valid/(length(xs)*length(ys)),
    mean_pointwise_pmc=mean(p),
    pmc_q025=unname(quantile(p, 0.025)), pmc_q975=unname(quantile(p, 0.975)),
    determinant_only_extra=sum(determinant > 0 & !valid))
}

estimate_one <- function(bounds, partial_r, rxz2, ryz2, nu, direction) {
  xs <- sort(runif(n_marginal, bounds[1], bounds[2]))
  ys <- sort(runif(n_marginal, bounds[1], bounds[2]))
  evaluate_pairs(xs, ys, partial_r, rxz2, ryz2, nu, direction)
}

# Reproduce the four original figure-script estimates BEFORE resetting the RNG.
# These single-run results are retained for comparison, not used as final labels.
original_single_run <- do.call(rbind, lapply(scenarios, function(s) {
  a <- estimate_one(s$bounds, s$partial_r, s$rxz2, s$ryz2, s$nu, s$direction)
  data.frame(scenario=s$name, as.list(a), check.names=FALSE)
}))
cat("\nOriginal final_ver single-run checks (original RNG order):\n")
print(original_single_run, digits=10, row.names=FALSE)

# Independent seed reset preserves every repetition of PMC_100run_diagnostics.R.
set.seed(123)
all_runs <- vector("list", length(scenarios))
summary_rows <- vector("list", length(scenarios))
for (i in seq_along(scenarios)) {
  s <- scenarios[[i]]
  rep_results <- replicate(n_repeats,
    estimate_one(s$bounds, s$partial_r, s$rxz2, s$ryz2, s$nu, s$direction))
  h <- as.numeric(rep_results["p_hat", ])
  n_accepted <- as.numeric(rep_results["n_admissible", ])
  all_runs[[i]] <- data.frame(
    scenario=s$name, run=seq_len(n_repeats), p_hat=h,
    n_admissible=n_accepted,
    acceptance_rate=as.numeric(rep_results["acceptance_rate", ]))
  summary_rows[[i]] <- data.frame(
    scenario=s$name,
    mean_accepted=mean(n_accepted), min_accepted=min(n_accepted),
    max_accepted=max(n_accepted), mean_p=mean(h), sd_across_runs=sd(h),
    se_mean_100=sd(h)/sqrt(n_repeats),
    first_50_mean=mean(h[1:50]), last_50_mean=mean(h[51:100]),
    half_difference=abs(mean(h[1:50])-mean(h[51:100])),
    n_repeats=n_repeats, n_marginal=n_marginal,
    candidate_pairs_per_run=n_marginal^2, threshold=threshold,
    partial_r=s$partial_r, rxz2=s$rxz2, ryz2=s$ryz2, nu=s$nu,
    direction=s$direction, mean_acceptance_rate=mean(n_accepted)/n_marginal^2,
    original_single_run_p=original_single_run$p_hat[i],
    figure_label_2dp=sprintf("%.2f", mean(h)),
    determinant_only_extra_pairs=sum(rep_results["determinant_only_extra", ]))
  cat("Completed", i, "/", length(scenarios), "scenarios:", s$name, "\n")
}
run_df <- do.call(rbind, all_runs)  # In memory for audit; no extra raw CSV output.
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(output_dir, "PMC_MC_100run_summary.csv"), row.names=FALSE)
cat("\nTable A3: conditional Pr(p_mc > 0.8 | admissible); mean of 100 run proportions.\n")
print(summary_df[, 1:10], digits=10, row.names=FALSE)
cat("\nSD is single-run simulation variability; MCSE of the 100-run mean = SD/10.\n",
    "The 1,000,000 Cartesian pairs per run are NOT independent observations.\n",
    "Admissible-pair counts are not effective sample sizes.\n",
    "First/last-half comparisons describe stability, not a convergence proof.\n",
    "Figure labels below use these 100-run means, rounded to two decimals.\n", sep="")
figure_probabilities <- summary_df$mean_p
retention_full_prob <- figure_probabilities[1L]
retention_restricted_prob <- figure_probabilities[2L]
education_full_prob <- figure_probabilities[3L]
education_restricted_prob <- figure_probabilities[4L]

make_surface <- function(lb, ub, rxy_z, rxz2, ryz2, nu, sign_r, n_grid = 801L) {
  x <- seq(lb, ub, length.out = n_grid)
  y <- seq(lb, ub, length.out = n_grid)
  X <- matrix(rep(x, times = n_grid), nrow = n_grid, ncol = n_grid)
  Y <- matrix(rep(y, each = n_grid), nrow = n_grid, ncol = n_grid)

  numerator <- rxy_z * sqrt((1 - ryz2) * (1 - rxz2)) - X * Y
  v_x <- 1 - rxz2 - X^2
  v_y <- 1 - ryz2 - Y^2
  determinant <- v_x * v_y - numerator^2

  # Conditions in manuscript Eq. (5): positive residual variances and determinant.
  admissible <- (v_x > 0) & (v_y > 0) & (determinant > 0)
  z <- matrix(NA_real_, nrow = n_grid, ncol = n_grid)
  z[admissible] <- pnorm(sign_r * numerator[admissible] /
                        sqrt(determinant[admissible] / nu) - 1.96)
  list(x = x, y = y, z = z, admissible = admissible)
}

plot_panel_bw <- function(surface, prior_title, reported_p, restricted_box = FALSE) {
  x <- surface$x
  y <- surface$y
  z <- surface$z
  admissible <- surface$admissible
  lo <- min(x)
  hi <- max(x)

  # Do not label contours on the curves; the separate legend is clearer.
  par(mar = c(3.8, 4.3, 3.5, 1.1), mgp = c(2.3, 0.7, 0), pty = 's')
  plot(NA_real_, NA_real_, type = 'n', xlim = c(lo, hi), ylim = c(lo, hi),
       xaxs = 'i', yaxs = 'i', asp = 1, axes = FALSE, xlab = '', ylab = '')

  # White feasible region / grey inadmissible region.
  image(x, y, 1 * admissible, add = TRUE, useRaster = TRUE,
        breaks = c(-0.5, 0.5, 1.5), col = c('grey88', 'white'))

  # Only two pmc levels: dashed 0.1, solid 0.8.
  z_range <- range(z, finite = TRUE)
  if (z_range[1] < 0.1 && z_range[2] > 0.1) {
    contour(x, y, z, add = TRUE, levels = 0.1,
            drawlabels = FALSE, col = 'grey45', lty = 2, lwd = 1.8)
  }
  if (z_range[1] < 0.8 && z_range[2] > 0.8) {
    contour(x, y, z, add = TRUE, levels = 0.8,
            drawlabels = FALSE, col = 'black', lty = 1, lwd = 2.3)
  }

  # Boundary of admissible region; not present if every point is admissible.
  if (any(admissible) && any(!admissible)) {
    contour(x, y, 1 * admissible, add = TRUE, levels = 0.5,
            drawlabels = FALSE, col = 'black', lty = 1, lwd = 1.3)
  }

  # Restricted range illustrated only on the full-range panel of retention.
  if (restricted_box) {
    rect(-0.5, -0.5, 0.5, 0.5, border = 'grey45', lty = 2, lwd = 1.2)
  }

  ticks <- if (lo == -1) seq(-1, 1, 0.5) else seq(-0.5, 0.5, 0.25)
  axis(1, at = ticks, cex.axis = 1.07)
  axis(2, at = ticks, las = 1, cex.axis = 1.07)
  box(lwd = 1)
  title(xlab = expression(r[xc]), ylab = expression(r[yc]),
        cex.lab = 1.25, line = 2.6)
  title(main = paste0(prior_title,
                      '\nPr(p_mc > 0.8) = ', sprintf('%.2f', reported_p)),
        cex.main = 1.12, font.main = 1, line = 0.6)
}

plot_legend_bw <- function(show_restricted_box) {
  par(mar = c(0, 0, 0, 0), pty = 'm')
  plot.new()
  entries <- c('p_mc = 0.8', 'p_mc = 0.1', 'Admissible boundary')
  line_types <- c(1, 2, 1)
  line_colors <- c('black', 'grey45', 'black')
  line_widths <- c(2.3, 1.8, 1.3)
  symbols <- c(NA, NA, NA)
  if (show_restricted_box) {
    entries <- c(entries, 'Restricted prior box')
    line_types <- c(line_types, 2)
    line_colors <- c(line_colors, 'grey45')
    line_widths <- c(line_widths, 1.2)
    symbols <- c(symbols, NA)
  }
  entries <- c(entries, 'Inadmissible region')
  line_types <- c(line_types, NA)
  line_colors <- c(line_colors, 'grey88')
  line_widths <- c(line_widths, NA)
  symbols <- c(symbols, 15)

  legend('center', legend = entries, ncol = 3, bty = 'n',
         lty = line_types, lwd = line_widths, col = line_colors,
         pch = symbols, pt.cex = 1.3,
         cex = 1.04, x.intersp = 0.8, y.intersp = 1.25)
}

plot_comparison_bw <- function(example, title_text, filename_stem,
                               baseline_prob, restricted_prob,
                               rxy_z, rxz2, ryz2, nu, sign_r,
                               show_restricted_box = FALSE) {
  # Precompute the full prior and the centre zoom for a direct comparison.
  full <- make_surface(-1, 1, rxy_z, rxz2, ryz2, nu, sign_r)
  restricted <- make_surface(-0.5, 0.5, rxy_z, rxz2, ryz2, nu, sign_r)

  draw <- function() {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE),
           heights = c(5.6, 1.0))
    par(oma = c(0.3, 0.2, 2.0, 0.2))

    plot_panel_bw(full, '(a) Prior: U(-1, 1)', baseline_prob,
                  restricted_box = show_restricted_box)
    plot_panel_bw(restricted, '(b) Prior: U(-0.5, 0.5)', restricted_prob)
    plot_legend_bw(show_restricted_box)

    mtext(title_text, side = 3, outer = TRUE, line = 0.5,
          cex = if (example == 1L) 1.45 else 1.22, font = 2)
  }

  # PDF: vector contours, labels, boundaries (scalable for LaTeX).
  pdf(file.path(output_dir, paste0(filename_stem, '.pdf')), width = 12, height = 7,
      family = 'Helvetica', useDingbats = FALSE)
  draw()
  dev.off()

  # PNG: 600 DPI, suitable for journal preview/submission.
  png(file.path(output_dir, paste0(filename_stem, '.png')), width = 12, height = 7,
      units = 'in', res = 600)
  draw()
  dev.off()

  if (interactive()) draw()
}
plot_comparison_bw(
  example = 1L,
  title_text = 'Kindergarten retention',
  filename_stem = 'Figure1_Retention_Priors_Comparison_BW_Legend',
  baseline_prob = retention_full_prob, restricted_prob = retention_restricted_prob,
  rxy_z = -0.13, rxz2 = 0.52^2, ryz2 = 0.45^2,
  nu = 7635, sign_r = -1,
  show_restricted_box = TRUE
)

plot_comparison_bw(
  example = 2L,
  title_text = "Father's occupation and educational attainment",
  filename_stem = 'Figure2_Education_Priors_Comparison_BW_Legend',
  baseline_prob = education_full_prob, restricted_prob = education_restricted_prob,
  # IMPORTANT: 0.224 and 0.177 are SQUARED multiple correlations.
  rxy_z = 0.29, rxz2 = 0.224, ryz2 = 0.177,
  nu = 10562, sign_r = +1,
  show_restricted_box = TRUE
)

cat('\nPDF and 600-DPI PNG figures saved to:', normalizePath(output_dir), '\n')

