# Run from repository root: Rscript tests/education.R [run1 run2 baseline_run]
# Full runs save list(runs=run_df, reconstruction=output, rng=.Random.seed) as audit.rds.
options(warn=2)
close <- function(a,b,tol=1e-12) {
  stopifnot(isTRUE(all.equal(a,b,tolerance=tol,check.attributes=FALSE)))
}
expressions <- as.list(parse('robustness.R'))
load_definition <- function(name) {
  definitions <- Filter(function(e) is.call(e) && identical(e[[1]],as.name('<-')) &&
                          identical(e[[2]],as.name(name)),expressions)
  stopifnot(length(definitions)==1L)
  eval(definitions[[1]],envir=.GlobalEnv)
}
# This assertion also fails on the restored pre-fix script's hardcoded .29.
if (any(vapply(expressions,function(e) is.call(e) &&
    identical(e[[1]],as.name('<-')) && identical(e[[2]],as.name('education_inputs')),logical(1)))) {
  load_definition('education_inputs')
}
load_definition('scenarios')
close(scenarios[[3]]$partial_r,.2984613894645406)
close(scenarios[[4]]$partial_r,.2984613894645406)
close(education_inputs$cross_term,.18838965505016483)
close(education_inputs$rxz2,.22352057677494624)
close(education_inputs$ryz2,.17686271729239786)
# Independently specified full Table 1 matrix and precision-matrix formula.
R <- matrix(c(1,.427,-.437,-.289,.427,1,-.317,-.350,
              -.437,-.317,1,.263,-.289,-.350,.263,1),4,byrow=TRUE)
stopifnot(min(eigen(R,symmetric=TRUE,only.values=TRUE)$values)>0)
K <- solve(R)
r <- -K[1,2]/sqrt(K[1,1]*K[2,2])
close(education_inputs$partial_r,r)
close(r*sqrt(10562/(1-r*r)),32.1381541016481)
stopifnot(scenarios[[3]]$nu==10562,scenarios[[4]]$nu==10562)
# Protected files must stay byte-identical to the restored starting commit.
protected <- c('retention.pdf','retention.png','minimum_effect.R',
               'minimum_effect_runs.csv','minimum_effect_summary.csv',
               'tests/minimum_effect.R','tests/minimum_effect_results.R')
stopifnot(system2('git',c('diff','--exit-code',
  '1cb146626f6a2ca057bbba8df8d44fdddc4a1b90','--',protected))==0)
cat('PASS: Education matrix, positive definiteness, precision partial, T0; protected files unchanged.\n')
args <- commandArgs(TRUE)
if (!length(args)) quit(status=0)
stopifnot(length(args)==3L)
s <- read.csv(file.path(args[1],'robustness_summary.csv'))
old <- read.csv(file.path(args[3],'robustness_summary.csv'))
for (name in c('robustness_summary.csv','education.png','retention.png')) {
  stopifnot(unname(tools::md5sum(file.path(args[1],name)))==
            unname(tools::md5sum(file.path(args[2],name))))
}
a <- readRDS(file.path(args[1],'audit.rds'))
b <- readRDS(file.path(args[2],'audit.rds'))
original <- readRDS(file.path(args[3],'audit.rds'))
stopifnot(identical(a,b),identical(a$rng,original$rng),
          identical(a$reconstruction,original$reconstruction),
          identical(a$runs[1:200,],original$runs[1:200,]),
          identical(s[1:2,],old[1:2,]),
          identical(readLines(file.path(args[1],'robustness_summary.csv'))[1:3],
                    readLines(file.path(args[3],'robustness_summary.csv'))[1:3]),
          unname(tools::md5sum(file.path(args[1],'retention.png')))==
          unname(tools::md5sum(file.path(args[3],'retention.png'))))
# Independently replay the original seed flow: two Retention priors consume
# the first 400,000 uniforms, followed by the two Education priors.
RNGkind('Mersenne-Twister','Inversion','Rejection')
set.seed(123)
invisible(runif(400000L))
V <- R[1:2,1:2]-R[1:2,3:4]%*%solve(R[3:4,3:4],R[3:4,1:2])
cutoff <- qnorm(.8)+1.96
for (j in 3:4) {
  bounds <- if(j==3) c(-1,1) else c(-.5,.5)
  for (i in 1:100) {
    x <- sort(runif(1000,bounds[1],bounds[2]))
    y <- sort(runif(1000,bounds[1],bounds[2]))
    x <- rep(x,times=1000); y <- rep(y,each=1000)
    vx <- V[1,1]-x*x; vy <- V[2,2]-y*y
    A <- V[1,2]-x*y; D <- vx*vy-A*A
    valid <- vx>0 & vy>0 & D>0
    tt <- A[valid]/sqrt(D[valid]/10562)
    stopifnot(all(is.finite(tt)),all(is.finite(pnorm(tt-1.96))))
    row <- a$runs[(j-1)*100+i,]
    close(row$n_admissible,sum(valid))
    close(row$p_hat,mean(tt>cutoff))
  }
  d <- a$runs[a$runs$scenario==s$scenario[j],]
  close(s$mean_p[j],mean(d$p_hat)); close(s$sd_across_runs[j],sd(d$p_hat))
  close(s$se_mean_100[j],sd(d$p_hat)/10)
  close(s$mean_accepted[j],mean(d$n_admissible))
  close(s$min_accepted[j],min(d$n_admissible))
  close(s$max_accepted[j],max(d$n_admissible))
}
cat('PASS: two full runs identical; all 200 Education repetitions independently verified; Retention runs, summary, PNG, reconstruction and RNG unchanged.\n')
