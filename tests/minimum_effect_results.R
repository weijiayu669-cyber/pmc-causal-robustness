# Run from repository root after two complete runs:
# Rscript tests/minimum_effect_results.R <run1_directory> <run2_directory>
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
options(warn=2)
for (name in c('minimum_effect_runs.csv','minimum_effect_summary.csv')) {
  paths <- file.path(args,name)
  stopifnot(unname(tools::md5sum(paths)[1])==unname(tools::md5sum(paths)[2]))
}
runs <- read.csv(file.path(args[1],'minimum_effect_runs.csv'))
s <- read.csv(file.path(args[1],'minimum_effect_summary.csv'))
stopifnot(nrow(runs)==200L,nrow(s)==2L,
          all(vapply(runs[vapply(runs,is.numeric,logical(1))],
                     function(x) all(is.finite(x)),logical(1))),
          all(runs$n_nonfinite==0),all(runs$min_D>0),all(runs$min_se>0),
          all(runs$max_relative_statistic_error<1e-12))
close <- function(x,y) stopifnot(isTRUE(all.equal(x,y,tolerance=1e-12,
                                               check.attributes=FALSE)))
for (i in seq_len(nrow(s))) {
  d <- runs[runs$prior==s$prior[i],]
  stopifnot(identical(d$run,1:100),all(d$n_candidate==1e6))
  close(d$Q0,d$n_above_Q0/d$n_admissible)
  close(d$Q_delta,d$n_above_Q_delta/d$n_admissible)
  close(d$difference_pp,100*(d$Q0-d$Q_delta))
  close(s$mean_Q0[i],mean(d$Q0)); close(s$mean_Q_delta[i],mean(d$Q_delta))
  close(s$sd_Q_delta[i],sd(d$Q_delta)); close(s$MCSE_Q_delta[i],sd(d$Q_delta)/10)
  close(s$mean_admissible[i],mean(d$n_admissible))
  close(s$difference_pp[i],100*(mean(d$Q0)-mean(d$Q_delta)))
  close(s$first_50_Q_delta[i],mean(d$Q_delta[1:50]))
  close(s$last_50_Q_delta[i],mean(d$Q_delta[51:100]))
}
# Compare full-precision baseline and diagnostics with the unchanged A.5 file.
a5 <- read.csv('robustness_summary.csv')[1:2,]
close(s$mean_Q0,a5$mean_p); close(s$sd_Q0,a5$sd_across_runs)
close(s$MCSE_Q0,a5$se_mean_100); close(s$mean_admissible,a5$mean_accepted)
close(s$first_50_Q0,a5$first_50_mean); close(s$last_50_Q0,a5$last_50_mean)
stopifnot(all(abs(s$mean_Q0-c(.678886,.836594))<5e-7))

# Replay the unchanged main evaluator for EVERY run using the same RNG stream.
exprs <- parse('robustness.R')
eval_expr <- Filter(function(e) is.call(e) && identical(e[[1]],as.name('<-')) &&
                     identical(e[[2]],as.name('evaluate_pairs')),as.list(exprs))
stopifnot(length(eval_expr)==1L)
eval(eval_expr[[1]])
threshold <- .8
RNGkind('Mersenne-Twister','Inversion','Rejection')
set.seed(123)
priors <- list(c(-1,1),c(-.5,.5))
for (j in seq_along(priors)) for (r in 1:100) {
  xs <- sort(runif(1000,priors[[j]][1],priors[[j]][2]))
  ys <- sort(runif(1000,priors[[j]][1],priors[[j]][2]))
  result <- evaluate_pairs(xs,ys,-.13,.52^2,(-.45)^2,7635,-1)
  row <- runs[(j-1)*100+r,]
  close(row$Q0,result['p_hat'])
  close(row$n_admissible,result['n_admissible'])
  if(r%%50==0) cat('Baseline replay prior',j,'run',r,'passed\n')
}
cat('PASS: byte-identical repeat CSVs, all finite diagnostics, summary arithmetic, A.5 statistics and all 200 baseline runs.\n')
cat('Maximum A10 absolute error:',max(runs$max_abs_statistic_error),'\n')
cat('Maximum A10 scaled error:',max(runs$max_relative_statistic_error),'\n')
print(s[,c('prior','first_50_Q_delta','last_50_Q_delta')],digits=15,row.names=FALSE)
cat('\n'); print(sessionInfo())
