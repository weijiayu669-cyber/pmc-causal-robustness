# Run from repository root after two complete runs:
# Rscript tests/robustness_results.R <run1_directory> <run2_directory>
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
options(warn=2)
close <- function(x,y,tol=1e-12) {
  stopifnot(isTRUE(all.equal(x,y,tolerance=tol,check.attributes=FALSE)))
}
for (name in c('robustness_runs.csv','robustness_summary.csv')) {
  paths <- file.path(args,name)
  stopifnot(all(file.exists(paths)),
            unname(tools::md5sum(paths)[1])==unname(tools::md5sum(paths)[2]))
}
runs <- read.csv(file.path(args[1],'robustness_runs.csv'))
s <- read.csv(file.path(args[1],'robustness_summary.csv'))
stopifnot(nrow(runs)==400L,nrow(s)==4L,
          all(vapply(runs[vapply(runs,is.numeric,logical(1))],
                     function(x) all(is.finite(x)),logical(1))),
          all(vapply(s[vapply(s,is.numeric,logical(1))],
                     function(x) all(is.finite(x)),logical(1))),
          all(runs$n_nonfinite==0),all(runs$min_D>0),
          all(runs$min_vx>0),all(runs$min_vy>0),
          all(runs$n_admissible>0),all(runs$n_admissible<=1e6),
          all(runs$n_above>=0),all(runs$n_above<=runs$n_admissible),
          all(runs$n_admissible==floor(runs$n_admissible)),
          all(runs$n_above==floor(runs$n_above)),
          all(s$n_repeats==100),all(s$n_marginal==1000),
          all(s$candidate_pairs_per_run==1e6),all(s$threshold==.8),
          all(s$determinant_only_extra_pairs==0))
expected_names <- c('Retention U(-1,1)','Retention U(-0.5,0.5)',
                    'Education U(-1,1)','Education U(-0.5,0.5)')
stopifnot(identical(s$scenario,expected_names),
          identical(runs$scenario,rep(expected_names,each=100)))
close(s$partial_r,c(-.13,-.13,rep(.2984613894645406,2)))
close(s$rxz2,c(.2704,.2704,rep(.22352057677494624,2)))
close(s$ryz2,c(.2025,.2025,rep(.17686271729239786,2)))
close(s$nu,c(7635,7635,10562,10562))
close(s$direction,c(-1,-1,1,1))
for (i in 1:4) {
  d <- runs[runs$scenario==s$scenario[i],]
  stopifnot(identical(d$run,1:100))
  close(d$p_hat,d$n_above/d$n_admissible)
  close(d$acceptance_rate,d$n_admissible/1e6)
  close(s$mean_p[i],mean(d$p_hat))
  close(s$sd_across_runs[i],sd(d$p_hat))
  close(s$se_mean_100[i],sd(d$p_hat)/10)
  close(s$first_50_mean[i],mean(d$p_hat[1:50]))
  close(s$last_50_mean[i],mean(d$p_hat[51:100]))
  close(s$half_difference[i],abs(mean(d$p_hat[1:50])-mean(d$p_hat[51:100])))
  close(s$mean_accepted[i],mean(d$n_admissible))
  close(s$min_accepted[i],min(d$n_admissible))
  close(s$max_accepted[i],max(d$n_admissible))
  close(s$mean_acceptance_rate[i],mean(d$acceptance_rate))
  stopifnot(sprintf('%.3f',s$figure_label_3dp[i])==sprintf('%.3f',mean(d$p_hat)))
}
# Historical Retention baseline stays unchanged after the Education correction.
stopifnot(all(abs(s$mean_p[1:2]-c(.678886,.836594))<5e-7))

# Independent seed replay: constants and threshold event do not load production
# inputs.R, evaluate_pairs(), or its pnorm-based classification.
# Education residual covariance comes directly from Table 1's control block.
RZZ <- matrix(c(1,.263,.263,1),2)
RxyZ <- rbind(c(-.437,-.289),c(-.317,-.350))
V_edu <- matrix(c(1,.427,.427,1),2) - RxyZ %*% solve(RZZ,t(RxyZ))
V_ret <- matrix(c(1-.2704,-.13*sqrt((1-.2704)*(1-.2025)),
                  -.13*sqrt((1-.2704)*(1-.2025)),1-.2025),2)
RNGkind('Mersenne-Twister','Inversion','Rejection')
set.seed(123)
cutoff <- qnorm(.8)+1.96
for (j in 1:4) {
  bounds <- if(j%%2==1) c(-1,1) else c(-.5,.5)
  V <- if(j<=2) V_ret else V_edu
  nu <- if(j<=2) 7635 else 10562
  direction <- if(j<=2) -1 else 1
  for (r in 1:100) {
    xs <- sort(runif(1000,bounds[1],bounds[2]))
    ys <- sort(runif(1000,bounds[1],bounds[2]))
    a <- rep(xs,times=1000); b <- rep(ys,each=1000)
    vx <- V[1,1]-a*a; vy <- V[2,2]-b*b
    cross <- V[1,2]-a*b
    D <- vx*vy-cross*cross
    valid <- vx>0 & vy>0 & D>0
    # Every admissible pair must have finite signed T and a finite probability.
    tt <- cross[valid]/sqrt(D[valid]/nu)
    p <- pnorm(direction*tt-1.96)
    stopifnot(all(is.finite(tt)),all(is.finite(p)),all(p>=0 & p<=1))
    above <- direction*cross[valid] > cutoff*sqrt(D[valid]/nu)
    row <- runs[(j-1)*100+r,]
    close(row$n_admissible,sum(valid))
    close(row$n_above,sum(above))
    close(row$p_hat,mean(above))
    # Independent matrix multiplication differs by roundoff near D=0.
    stopifnot(abs(row$min_D-min(D[valid]))<1e-14)
    close(row$min_vx,min(vx[valid]))
    close(row$min_vy,min(vy[valid]))
    if(r%%50==0) cat('Independent scenario',j,'run',r,'passed\n')
  }
}
cat('PASS: byte-identical CSVs; 400 independent seed replays; all admissible pairs finite; corrected Education, retained Retention, counts, means, SD and MCSE=SD/10.\n')
print(s[,c('scenario','mean_p','sd_across_runs','se_mean_100',
            'first_50_mean','last_50_mean')],digits=15,row.names=FALSE)
cat('\n'); print(sessionInfo())
