# Appendix A.6: Hong and Raudenbush (2005) retention illustration ONLY.
# Base R; Rscript minimum_effect.R [output_directory], or source().
# Extends the supplied PMC working OLS model; does NOT refit the original HLM.
# Original correlations, strict admissibility, df, cutoff, draws and RNG retained.
# Main Delta uses the approximate literature conversion, not an exact known SD.
# Table A4 is written to minimum_effect_summary.csv; A6_summary.csv
# additionally retains all secondary criteria and the working-SD scale check.

working_scale <- function() {
  # Verbatim size/meany/sdy inputs from the supplied robustness.R mcdat().
  # These are working summaries; original reading missingness is not corrected.
  cells <- data.frame(
    n=c(9,12,14,12,24,23,47,48,46,48,47,49,45,39,8,
        3054,1661,983,323,441,154,213,143,85,49,35,16,8,3),
    mean=c(47.41,47.25,42.66,45.24,44.07,36.63,37.41,39.79,37.78,
           35,34.3,32.06,31.09,33.87,32.84,59.22,55.62,52.49,
           52.45,50.05,48.77,45.89,46,43.28,42.18,40.45,33.03,35.27,40.47),
    sd=c(11.73,15.50,15.92,17.32,14.43,13.88,10.92,15.37,12.27,
         9.77,9.05,9.31,11.96,9.15,7.71,11.64,11.30,11.95,12.53,
         12.15,12.67,12.40,13.37,11.52,13.34,13.08,14.33,8.25,8.60))
  n <- sum(cells$n)
  ybar <- sum(cells$n*cells$mean)/n
  sy <- sqrt(sum((cells$n-1)*cells$sd^2 + cells$n*(cells$mean-ybar)^2)/(n-1))
  sx <- sqrt(n/(n-1)*(471/n)*(7168/n)) # binary X: 1=retained, 0=promoted
  list(n=n, n_retained=471, n_promoted=7168, mean_y=ybar,
       sd_y=sy, sd_x=sx, scale_ratio=sy/sx)
}

pointwise <- function(a, b, scale_ratio, delta) {
  # C is orthogonal to the observed Z, as specified in Appendix A.2.
  q <- -0.13*sqrt((1-0.52^2)*(1-0.45^2))
  vx <- 1-0.52^2-a*a
  vy <- 1-0.45^2-b*b
  A <- q-a*b
  D <- vx*vy-A*A
  valid <- vx>0 & vy>0 & D>0
  beta <- se <- tt <- p0 <- pd <- pj <- rep(NA_real_,length(A))
  beta[valid] <- scale_ratio*A[valid]/vx[valid]
  se[valid] <- scale_ratio*sqrt(D[valid]/7635)/vx[valid]
  tt[valid] <- A[valid]/sqrt(D[valid]/7635)
  # Original formula, unchanged. It is directional normal-approximation power.
  p0[valid] <- pnorm(-tt[valid]-1.96)
  # H0: beta >= -Delta vs H1: beta < -Delta, same 1.96 cutoff (alpha~.025).
  pd[valid] <- pnorm((-delta-beta[valid])/se[valid]-1.96)
  # Separate sampling event: B < -1.96*SE AND B < -Delta,
  # where B ~ N(beta, SE^2), treating the scenario-specific SE as fixed.
  pj[valid] <- pnorm(-beta[valid]/se[valid]-pmax(1.96,delta/se[valid]))
  list(valid=valid,beta=beta,se=se,t=tt,p_zero=p0,p_delta=pd,p_joint=pj)
}

output_dir <- local({
  files <- Filter(function(x) is.character(x) && length(x)==1,
                  lapply(sys.frames(),function(x) x$ofile))
  arg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
  self <- if(length(files)) tail(files,1)[[1]] else if(length(arg)) sub('^--file=','',arg[1]) else NULL
  extra <- if(length(files)) character() else commandArgs(TRUE)
  if(length(extra)) extra[1] else if(is.null(self)) getwd() else dirname(normalizePath(self))
})
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
sc <- working_scale()
delta <- 0.20*9.01/0.67 # 2.689552..., approximate literature-referenced points
delta_working_sd <- 0.20*sc$sd_y # secondary check, explicitly a different SD basis
write.csv(data.frame(as.list(sc),literature_reference_sd=9.01/.67,
                     delta_main_points=delta,delta_working_sd_points=delta_working_sd,
                     main_delta_in_working_sd=delta/sc$sd_y),
          file.path(output_dir,'scale_audit.csv'),row.names=FALSE)

RNGkind('Mersenne-Twister','Inversion','Rejection')
set.seed(123)
# The diagnostics stream in robustness.R resets seed=123 before the two
# retention scenarios. Its earlier mcdat() reconstruction consumes no draws
# from this reset stream; omitting it here does not alter the A.5 results.
n_marginal <- 1000L
n_runs <- 100L
priors <- list('U(-1,1)'=c(-1,1),'U(-0.5,0.5)'=c(-.5,.5))
criteria <- c('original_NHST_power','minimum_effect_directional_power',
              'joint_sampling_event_power','NHST_power_and_plugin_effect_filter',
              'directional_power_at_0.20_working_SD')
rows <- vector('list',length(priors)*n_runs*length(criteria))
axes <- vector('list',length(priors)*n_runs)
k <- 0L
for (j in seq_along(priors)) {
  for (r in seq_len(n_runs)) {
    xs <- sort(runif(n_marginal,priors[[j]][1],priors[[j]][2]))
    ys <- sort(runif(n_marginal,priors[[j]][1],priors[[j]][2]))
    axes[[(j-1L)*n_runs+r]] <- list(prior=names(priors)[j],run=r,xs=xs,ys=ys)
    a <- rep(xs,times=n_marginal)
    b <- rep(ys,each=n_marginal)
    z <- pointwise(a,b,sc$scale_ratio,delta)
    v <- z$valid
    p0 <- z$p_zero[v]; pd <- z$p_delta[v]; pj <- z$p_joint[v]
    beta <- z$beta[v]; se <- z$se[v]
    pw <- pnorm((-delta_working_sd-beta)/se-1.96)
    # Identities / inequalities apply to every admissible configuration.
    stopifnot(max(abs(beta/se-z$t[v]))<1e-8,
              all(pd<=pj+1e-12), all(pj<=p0+1e-12),
              all(pd<=pw+1e-12))
    flags <- list(p0>.8,pd>.8,pj>.8,(p0>.8)&(-beta>=delta),pw>.8)
    powers <- list(p0,pd,pj,NA_real_,pw)
    for (h in seq_along(criteria)) {
      k <- k+1L
      rows[[k]] <- data.frame(prior=names(priors)[j],run=r,criterion=criteria[h],
        scenario_proportion=mean(flags[[h]]),n_above=sum(flags[[h]]),
        n_admissible=sum(v),n_candidate=n_marginal^2,
        acceptance_rate=mean(v),mean_pointwise_power=mean(powers[[h]]),
        delta_points=if(h==1) 0 else if(h==5) delta_working_sd else delta)
    }
    if(r%%10L==0L) cat(names(priors)[j],': completed',r,'of',n_runs,'runs\n')
  }
}
runs <- do.call(rbind,rows)
write.csv(runs,file.path(output_dir,'A6_runs_100.csv'),row.names=FALSE)
saveRDS(axes,file.path(output_dir,'A6_sampled_axes.rds'),compress='gzip')
summary_rows <- list()
for (p in names(priors)) for (c in criteria) {
  d <- runs[runs$prior==p & runs$criterion==c,]
  s <- sd(d$scenario_proportion)
  summary_rows[[length(summary_rows)+1L]] <- data.frame(
    prior=p,criterion=c,mean=mean(d$scenario_proportion),sd_across_runs=s,
    MCSE=s/sqrt(n_runs),mean_admissible=mean(d$n_admissible),
    sd_admissible=sd(d$n_admissible),min_admissible=min(d$n_admissible),
    max_admissible=max(d$n_admissible),total_admissible=sum(d$n_admissible),
    mean_acceptance_rate=mean(d$acceptance_rate),
    first_50_mean=mean(d$scenario_proportion[1:50]),
    last_50_mean=mean(d$scenario_proportion[51:100]),
    mean_pointwise_power=mean(d$mean_pointwise_power),delta_points=d$delta_points[1],
    n_runs=n_runs,n_marginal=n_marginal,candidate_pairs_per_run=n_marginal^2)
}
summary_df <- do.call(rbind,summary_rows)
write.csv(summary_df,file.path(output_dir,'A6_summary.csv'),row.names=FALSE)
write.csv(summary_df[summary_df$criterion %in% criteria[1:2],],
          file.path(output_dir,'minimum_effect_summary.csv'),row.names=FALSE)

old <- summary_df[summary_df$criterion==criteria[1],]
reference <- data.frame(prior=names(priors),A5_mean=c(.678886,.836594),
                        A5_SD=c(.008726,.006524),A5_MCSE=c(.000873,.000652),
                        A5_mean_admissible=c(594284.2,1000000),
                        A5_first50=c(.679303,.836014),A5_last50=c(.678469,.837174))
baseline <- data.frame(reference,reproduced_mean=old$mean,
    reproduced_SD=old$sd_across_runs,reproduced_MCSE=old$MCSE,
    reproduced_mean_admissible=old$mean_admissible,
    reproduced_first50=old$first_50_mean,reproduced_last50=old$last_50_mean)
baseline$mean_difference <- baseline$reproduced_mean-baseline$A5_mean
baseline$matches_A5_rounding <- abs(baseline$mean_difference)<5e-7 &
    abs(baseline$reproduced_SD-baseline$A5_SD)<5e-7 &
    abs(baseline$reproduced_MCSE-baseline$A5_MCSE)<5e-7 &
    abs(baseline$reproduced_mean_admissible-baseline$A5_mean_admissible)<.05 &
    abs(baseline$reproduced_first50-baseline$A5_first50)<5e-7 &
    abs(baseline$reproduced_last50-baseline$A5_last50)<5e-7
write.csv(baseline,file.path(output_dir,'A5_baseline_reproduction.csv'),row.names=FALSE)
stopifnot(all(baseline$matches_A5_rounding))

comparisons <- list()
for(p in names(priors)) for(c in criteria[-1]) {
  z <- runs[runs$prior==p & runs$criterion==c,]
  b <- runs[runs$prior==p & runs$criterion==criteria[1],]
  diff <- z$scenario_proportion-b$scenario_proportion
  comparisons[[length(comparisons)+1L]] <- data.frame(prior=p,criterion=c,
    original_mean=mean(b$scenario_proportion),new_mean=mean(z$scenario_proportion),
    difference_new_minus_old=mean(diff),paired_difference_SD=sd(diff),
    paired_difference_MCSE=sd(diff)/sqrt(n_runs))
}
write.csv(do.call(rbind,comparisons),file.path(output_dir,'A6_comparison.csv'),row.names=FALSE)
writeLines(c('Analysis: retention PMC working OLS model only.',
  'Approximate reference: Delta=0.20*(9.01/0.67) reading-score points.',
  'Main effect/SE scale: deterministic pooled sample moments of supplied summaries.',
  'MCSE=SD across 100 independent axis draws / sqrt(100); no iid-pair SE.',
  capture.output(sessionInfo())),file.path(output_dir,'sessionInfo.txt'))
cat('\nA.5 reproduction:\n'); print(baseline,row.names=FALSE,digits=10)
cat('\nA.6: scenario proportions conditional on admissibility:\n')
print(summary_df[,c('prior','criterion','mean','sd_across_runs','MCSE','mean_admissible')],
      row.names=FALSE,digits=10)
