# Run from repository root: Rscript tests/retention_audit.R
# Frozen pre-extraction implementation: exact draws AND RNG state are contractual.
vrnorm <- Vectorize(rnorm,c("n","mean","sd"))

legacy_mcdat <- function() {
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


source('retention_reconstruction.R')
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
for (seed in c(1L, 123L, 456L)) {
  set.seed(seed)
  expected <- replicate(3L, legacy_mcdat())
  expected_state <- .Random.seed
  set.seed(seed)
  observed <- replicate(3L, mcdat())
  stopifnot(identical(expected, observed), identical(expected_state, .Random.seed))
}
cells <- retention_cells()
stopifnot(nrow(cells)==29L, sum(cells$n)==7639L,
          sum(cells$n[cells$x==1])==471L, sum(cells$n[cells$x==0])==7168L,
          all(cells$sdz>0), all(cells$sdy>0))
set.seed(123)
independent <- mcdat(return_data=TRUE)
set.seed(123)
positive <- mcdat(.5, return_data=TRUE)
set.seed(123)
negative <- mcdat(-.5, return_data=TRUE)
# Common random numbers keep z unchanged, and correlate standardized errors.
k <- independent$strata
ze <- (independent$z-cells$meanz[k])/cells$sdz[k]
ye <- (independent$y-cells$meany[k])/cells$sdy[k]
stopifnot(identical(independent$z,positive$z), identical(independent$z,negative$z),
  max(abs(positive$y-(cells$meany[k]+cells$sdy[k]*(.5*ze+sqrt(.75)*ye))))<1e-12,
  max(abs(negative$y-(cells$meany[k]+cells$sdy[k]*(-.5*ze+sqrt(.75)*ye))))<1e-12)
stopifnot(abs(cor(ze,(positive$y-cells$meany[k])/cells$sdy[k])-.5)<.03,
          abs(cor(ze,(negative$y-cells$meany[k])/cells$sdy[k])+.5)<.03)
# Gaussian draws preserve distributional marginals; do not force sample moments.
stopifnot(abs(mean(ze))<.04, abs(sd(ze)-1)<.04,
          abs(mean(ye))<.04, abs(sd(ye)-1)<.04,
          mean(independent$z[k==1])!=cells$meanz[1])
matched <- retention_cells('outcome_matched')
stopifnot(sum(matched$n)==7608L, sum(cells$n-matched$n)==31L,
          all(matched$n<=cells$n),
          identical(cells[c('x','strata','meanz','sdz','meany','sdy')],
                    matched[c('x','strata','meanz','sdz','meany','sdy')]))
set.seed(123)
dat <- mcdat(sample='outcome_matched', return_data=TRUE)
stopifnot(nrow(dat)==7608L, all(is.finite(as.matrix(dat))),
          identical(as.integer(table(dat$strata)),as.integer(matched$n)))
for (bad in c(-1.1,1.1,NA,Inf)) {
  stopifnot(inherits(try(mcdat(bad),silent=TRUE),'try-error'))
}
cat('PASS: exact legacy draws/RNG; cell totals; missing outcomes; Gaussian marginals; within-cell sensitivity.\n')

# Independent transcription of the original publication (Tables 8/9, p.216).
source_cells <- read.csv('retention_source_cells.csv')
stopifnot(nrow(source_cells)==30L,source_cells$table8_n[30L]==0L,
          source_cells$table9_n[30L]==0L)
source_cells <- source_cells[source_cells$table8_n>0,]
stopifnot(all(cells$n==source_cells$table8_n),
          all(matched$n==source_cells$table9_n),
          all(cells$x==source_cells$x),
          all(cells$meanz==source_cells$table8_logit_mean),
          all(cells$sdz==source_cells$table8_logit_sd),
          all(cells$meany==source_cells$table9_reading_mean),
          all(cells$sdy==source_cells$table9_reading_sd),
          all(cells$n-matched$n==source_cells$count_deficit))
source('retention_audit.R')
audit_cells <- retention_cell_source_audit()
stopifnot(all(audit_cells$missing_reading==source_cells$count_deficit),
          sum(audit_cells$missing_reading)==31L)
# Partial-r/T identity independently checked through the direct OLS fit.
for (rho in c(-.5,0,.5)) {
  set.seed(123)
  dat <- mcdat(rho,return_data=TRUE)
  tt_lm <- unname(coef(summary(lm(y~x+z,data=dat)))['x','t value'])
  # OLS fit has N-3 df; the requested augmented PMC reference uses N-4.
  tt_requested <- tt_lm*sqrt((nrow(dat)-4)/(nrow(dat)-3))
  set.seed(123)
  rr <- mcdat(rho)[4L]
  stopifnot(abs(rr*sqrt((nrow(dat)-4)/(1-rr^2))-tt_requested)<1e-10)
}
# Check diagnostics against population moments and exact nonlinear averages.
# This small test is independent of the 1,000-run deliverable audit.
out_dir <- tempfile('retention-audit-')
dir.create(out_dir)
result <- run_retention_audit(out_dir,n_sim=100L)
for (key in c('rxz','ryz','rxy','partial_r')) {
  stopifnot(all(abs(result[[paste0('mean_',key)]][1:10]-
                    result[[paste0('analytic_',key)]][1:10]) <
                  6*result[[paste0('mcse_',key)]][1:10]))
}
set.seed(123)
rr <- replicate(100L,mcdat())[4L,]
tt <- rr*sqrt(7635/(1-rr^2))
ref0 <- result[result$scenario=='reference' & !is.na(result$rho_within) &
               result$rho_within==0,]
stopifnot(abs(ref0$mean_T0-mean(tt))<1e-13,
          abs(ref0$mean_abs_T0_discrepancy-mean(abs(tt+13.27)))<1e-13,
          abs(ref0$mcse_T0-sd(tt)/10)<1e-13,
          all(result$nu[1:5]==7635),all(result$nu[6:10]==7604),
          result$nu[11]==7635,result$mean_partial_r[11]==-.13)
unlink(out_dir,recursive=TRUE)
cat('PASS: source-table transcription; mixture oracle; OLS T identity; per-run nonlinear T and MCSE.\n')
