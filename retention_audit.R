# Standalone, assumption-dependent reconstruction audit; base R only.
# Rscript retention_audit.R (outputs beside this script); no formal PMC changes.
retention_audit_dir <- local({
  source_files <- lapply(sys.frames(), function(frame) frame$ofile)
  source_files <- Filter(function(path) is.character(path) && length(path)==1L,
                         source_files)
  script_arg <- grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)
  script_file <- if (length(source_files)) tail(source_files,1L)[[1L]] else
    if (length(script_arg)) sub('^--file=','',script_arg[1L]) else NULL
  if (is.null(script_file)) getwd() else dirname(normalizePath(script_file))
})
source(file.path(retention_audit_dir,'retention_reconstruction.R'))
source(file.path(retention_audit_dir,'inputs.R'))

# Population mixture moments provide an independent infinite-sample check.
# They do not force finite synthetic samples to equal the published moments.
retention_moments <- function(sample='reference',rho_within=0) {
  a <- retention_cells(sample)
  w <- a$n/sum(a$n)
  mx <- sum(w*a$x); mz <- sum(w*a$meanz); my <- sum(w*a$meany)
  vx <- sum(w*(a$x-mx)^2)
  vz <- sum(w*(a$sdz^2+(a$meanz-mz)^2))
  vy <- sum(w*(a$sdy^2+(a$meany-my)^2))
  rxz <- sum(w*(a$x-mx)*(a$meanz-mz))/sqrt(vx*vz)
  rxy <- sum(w*(a$x-mx)*(a$meany-my))/sqrt(vx*vy)
  ryz <- sum(w*(rho_within*a$sdz*a$sdy+
                (a$meanz-mz)*(a$meany-my)))/sqrt(vz*vy)
  c(rxz=rxz,ryz=ryz,rxy=rxy,
    partial_r=(rxy-rxz*ryz)/sqrt((1-rxz^2)*(1-ryz^2)))
}

retention_cell_source_audit <- function() {
  a <- retention_cells()
  b <- retention_cells('outcome_matched')
  data.frame(a[c('x','strata')],
    publication_stratum=c(seq_len(15L),seq_len(14L)),
    n_reference=a$n,n_outcome_matched=b$n,missing_reading=a$n-b$n,
    a[c('meanz','sdz','meany','sdy')],
    count_reference='Table 8, p. 216',count_outcome='Table 9, p. 216',
    z_summary_source='Table 8, p. 216',y_summary_source='Table 9, p. 216',
    matched_z_assumption='Full-cell logit moments reused; subset moments unknown',
    source_url=paste0('https://www.academia.edu/29274731/',
      'Effects_of_Kindergarten_Retention_Policy_on_Childrens_Cognitive_Growth_in_Reading_and_Mathematics'))
}

run_retention_audit <- function(output_dir=retention_audit_dir,n_sim=1000L,seed=123L) {
  stopifnot(length(n_sim)==1L,n_sim>=2L,n_sim==as.integer(n_sim))
  RNGkind('Mersenne-Twister','Inversion','Rejection')
  # Fixed sensitivity grid, never fitted or chosen to match the reported -13.27.
  rho_grid <- c(-.5,-.25,0,.25,.5)
  rows <- list()
  for (sample in c('reference','outcome_matched')) {
    n <- sum(retention_cells(sample)$n)
    nu <- n-4L
    for (rho in rho_grid) {
      # Reset per scenario; the same latent normals are used across rho values
      # within each count sample. Unlike robustness.R, no initial draw is skipped.
      set.seed(seed)
      sims <- t(replicate(n_sim,mcdat(rho,sample)))
      colnames(sims) <- c('rxz','ryz','rxy','partial_r')
      # Augmented PMC reference uses N-4 df. A direct reduced lm(y~x+z)
      # has N-3 df and therefore a different t statistic.
      tt <- sims[,'partial_r']*sqrt(nu/(1-sims[,'partial_r']^2))
      discrepancy <- abs(tt-(-13.27))
      stats <- cbind(sims,T0=tt,abs_T0_discrepancy=discrepancy)
      summary <- c(setNames(colMeans(stats),paste0('mean_',colnames(stats))),
        setNames(apply(stats,2L,sd),paste0('sd_',colnames(stats))),
        setNames(apply(stats,2L,sd)/sqrt(n_sim),paste0('mcse_',colnames(stats))))
      analytic <- retention_moments(sample,rho)
      rows[[length(rows)+1L]] <- data.frame(
        scenario=sample,rho_within=rho,N=n,nu=nu,n_sim=n_sim,seed=seed,
        source='Assumption-dependent Gaussian cell reconstruction',
        z_subset_assumption=if (sample=='outcome_matched')
          'Table 8 full-cell moments reused for outcome-observed subset' else
          'Table 8 counts with Table 9 outcome summaries',
        as.list(summary),
        abs_mean_T0_discrepancy=abs(mean(tt)-(-13.27)),
        analytic_rxz=analytic['rxz'],analytic_ryz=analytic['ryz'],
        analytic_rxy=analytic['rxy'],analytic_partial_r=analytic['partial_r'],
        analytic_T0=analytic['partial_r']*sqrt(nu/(1-analytic['partial_r']^2)),
        row.names=NULL)
    }
  }
  diagnostics <- do.call(rbind,rows)
  # Explicit comparison row, derived from the retained rounded working inputs.
  # These values are neither estimated nor selected by this sensitivity audit.
  formal <- diagnostics[1L,]
  formal[] <- NA
  formal$scenario <- 'formal_rounded_working_inputs'
  formal$N <- 7639L; formal$nu <- retention_inputs$nu; formal$n_sim <- 0L
  formal$source <- 'Fixed rounded PMC inputs; not reconstructed or recalibrated'
  formal$z_subset_assumption <- 'Not applicable'
  formal$mean_rxz <- sqrt(retention_inputs$rxz2)
  # Negative y-z sign is established by the reconstruction/source contrast.
  formal$mean_ryz <- -sqrt(retention_inputs$ryz2)
  formal$mean_partial_r <- retention_inputs$partial_r
  formal$mean_rxy <- retention_inputs$partial_r*
    sqrt((1-retention_inputs$rxz2)*(1-retention_inputs$ryz2))+
    formal$mean_rxz*formal$mean_ryz
  formal$mean_T0 <- retention_inputs$partial_r*
    sqrt(retention_inputs$nu/(1-retention_inputs$partial_r^2))
  formal$mean_abs_T0_discrepancy <- abs(formal$mean_T0-(-13.27))
  formal$abs_mean_T0_discrepancy <- formal$mean_abs_T0_discrepancy
  diagnostics <- rbind(diagnostics,formal)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  write.csv(diagnostics,file.path(output_dir,'retention_diagnostics.csv'),row.names=FALSE)
  write.csv(retention_cell_source_audit(),
            file.path(output_dir,'retention_cell_source_audit.csv'),row.names=FALSE)
  print(diagnostics[c('scenario','rho_within','N','nu','mean_partial_r',
                      'mean_T0','mcse_T0','abs_mean_T0_discrepancy')],row.names=FALSE)
  cat('Assumption-dependent diagnostics only; formal PMC inputs unchanged.\n')
  invisible(diagnostics)
}

if (sys.nframe()==0L) run_retention_audit()
