# Appendix A.6: Hong and Raudenbush (2005) retention illustration only.
# Base R; Rscript minimum_effect.R [output_directory], or source().
# A.6 uses the fixed numerical cutoff 0.20*13.53 in place of 1.96.
# The original NHST baseline, admissibility, draws and RNG are unchanged.
# Outputs: minimum_effect_runs.csv and minimum_effect_summary.csv.

pointwise <- function(a, b, d_fixed) {
  # C is orthogonal to observed Z, as in the existing working OLS model.
  q <- -0.13*sqrt((1-0.52^2)*(1-(-0.45)^2))
  vx <- 1-0.52^2-a*a
  vy <- 1-(-0.45)^2-b*b
  A <- q-a*b
  D <- vx*vy-A*A
  valid <- vx>0 & vy>0 & D>0
  tt <- p0 <- p_fixed <- rep(NA_real_,length(A))
  tt[valid] <- A[valid]/sqrt(D[valid]/7635)
  p0[valid] <- pnorm(-tt[valid]-1.96)
  p_fixed[valid] <- pnorm(-tt[valid]-d_fixed)
  list(valid=valid,t=tt,p_zero=p0,p_fixed=p_fixed,D=D)
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
# Supplied source: Hong & Raudenbush (2005), Table 2A, reading SD=13.53.
# Kraft (2020) motivates 0.20 SD; here its numerical value replaces d_alpha.
# This is a fixed-cutoff analysis, not a reading-score exceedance probability.
d_fixed <- 0.20 * 13.53
RNGkind('Mersenne-Twister','Inversion','Rejection')
set.seed(123)
# Match A.5: reset once before the full prior, then continue into restricted.
n_marginal <- 1000L
n_runs <- 100L
priors <- list('U(-1,1)'=c(-1,1),'U(-0.5,0.5)'=c(-.5,.5))
rows <- vector('list',length(priors)*n_runs)
for (j in seq_along(priors)) {
  for (r in seq_len(n_runs)) {
    xs <- sort(runif(n_marginal,priors[[j]][1],priors[[j]][2]))
    ys <- sort(runif(n_marginal,priors[[j]][1],priors[[j]][2]))
    a <- rep(xs,times=n_marginal)
    b <- rep(ys,each=n_marginal)
    z <- pointwise(a,b,d_fixed)
    v <- z$valid
    stopifnot(!anyNA(v),any(v))
    p0 <- z$p_zero[v]; p_fixed <- z$p_fixed[v]
    n_nonfinite <- sum(!is.finite(z$t[v]) | !is.finite(p0) | !is.finite(p_fixed))
    # Do not silently discard NA/NaN/Inf from admissible configurations.
    stopifnot(n_nonfinite==0L,all(p0>=0 & p0<=1),
              all(p_fixed>=0 & p_fixed<=1),all(p_fixed<=p0))
    Q0 <- mean(p0>.8)
    Q_fixed <- mean(p_fixed>.8)
    rows[[(j-1L)*n_runs+r]] <- data.frame(
      prior=names(priors)[j],run=r,Q0=Q0,Q_fixed=Q_fixed,
      difference_pp=100*(Q0-Q_fixed),n_above_Q0=sum(p0>.8),
      n_above_Q_fixed=sum(p_fixed>.8),n_admissible=sum(v),
      n_candidate=n_marginal^2,mean_pointwise_original=mean(p0),
      mean_pointwise_fixed=mean(p_fixed),d_fixed=d_fixed,min_D=min(z$D[v]),
      n_nonfinite=n_nonfinite)
    if(r%%10L==0L) cat(names(priors)[j],': completed',r,'of',n_runs,'runs\n')
  }
}
runs <- do.call(rbind,rows)
write.csv(runs,file.path(output_dir,'minimum_effect_runs.csv'),row.names=FALSE)
summary_rows <- lapply(names(priors),function(p) {
  d <- runs[runs$prior==p,]
  data.frame(prior=p,mean_Q0=mean(d$Q0),mean_Q_fixed=mean(d$Q_fixed),
    sd_Q_fixed=sd(d$Q_fixed),MCSE_Q_fixed=sd(d$Q_fixed)/sqrt(n_runs),
    mean_admissible=mean(d$n_admissible),difference_pp=mean(d$difference_pp),
    sd_Q0=sd(d$Q0),MCSE_Q0=sd(d$Q0)/sqrt(n_runs),
    sd_difference_pp=sd(d$difference_pp),
    MCSE_difference_pp=sd(d$difference_pp)/sqrt(n_runs),
    first_50_Q0=mean(d$Q0[1:50]),last_50_Q0=mean(d$Q0[51:100]),
    first_50_Q_fixed=mean(d$Q_fixed[1:50]),last_50_Q_fixed=mean(d$Q_fixed[51:100]),
    min_admissible=min(d$n_admissible),max_admissible=max(d$n_admissible),
    d_fixed=d_fixed,n_runs=n_runs,n_marginal=n_marginal,
    candidate_pairs_per_run=n_marginal^2)
})
summary_df <- do.call(rbind,summary_rows)
write.csv(summary_df,file.path(output_dir,'minimum_effect_summary.csv'),row.names=FALSE)
cat('\nFixed cutoff =',format(d_fixed,digits=17),'\n')
print(summary_df,digits=15,row.names=FALSE)
