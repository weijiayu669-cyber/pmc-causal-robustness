# Run from repository root: Rscript tests/minimum_effect.R
# Load the real evaluator and cutoff without running the Monte Carlo driver.
exprs <- as.list(parse('minimum_effect.R'))
for (name in c('pointwise','sx','sy','delta_min')) {
  definition <- Filter(function(e) is.call(e) && identical(e[[1]],as.name('<-')) &&
                        identical(e[[2]],as.name(name)),exprs)
  if(length(definition)!=1L) stop('Missing effect-size definition: ',name)
  eval(definition[[1]])
}
options(warn=2)
stopifnot(abs(delta_min-2.706)<1e-14,sx==.240548,sy==13.53)

# Independent T-statistic oracle: exact covariance and a direct OLS fit.
n <- 7639L
k <- seq_len(n)
basis <- qr.Q(qr(cbind(1,sin(k),cos(k),sin(2*k),cos(2*k))))[,2:5]*sqrt(n-1)
rxz <- .52; ryz <- -.45
rxy <- -.13*sqrt((1-rxz^2)*(1-ryz^2))+rxz*ryz
for (ab in list(c(0,0),c(.3,-.3),c(-.4,-.3))) {
  a <- ab[1]; b <- ab[2]
  R <- matrix(c(1,rxy,rxz,a,rxy,1,ryz,b,rxz,ryz,1,0,a,b,0,1),4)
  dat <- as.data.frame(basis %*% chol(R))
  names(dat) <- c('x','y','z','c')
  dat$x <- dat$x*sx; dat$y <- dat$y*sy
  fit <- lm(y~x+z+c,data=dat)
  tt <- coef(summary(fit))['x','t value']
  se <- coef(summary(fit))['x','Std. Error']
  beta <- unname(coef(fit)['x'])
  z <- pointwise(a,b,sy/sx,delta_min)
  stopifnot(z$valid,df.residual(fit)==7635,abs(z$t-tt)<1e-10,
            abs(z$se_f-se)<1e-10,abs(z$t-beta/z$se_f)<1e-10,
            abs(z$d_delta-2.706/se)<1e-10,
            abs(z$p_zero-pnorm(-tt-1.96))<1e-10,
            abs(z$p_delta-pnorm((-2.706-beta)/se))<1e-10)
  # A cutoff of 1.96*SE must recover the original NHST probability.
  original <- pointwise(a,b,sy/sx,1.96*se)
  stopifnot(identical(z$t,original$t),identical(z$valid,original$valid),
            abs(original$p_delta-original$p_zero)<1e-10)
  # At beta=-Delta, probability is 0.5: detects an extra -1.96.
  stopifnot(abs(pointwise(a,b,sy/sx,-beta)$p_delta-.5)<1e-10)
  # Rescaling scores and Delta together must not change the probability.
  scaled <- pointwise(a,b,2*sy/sx,2*delta_min)
  stopifnot(identical(z$p_delta,scaled$p_delta),
            identical(scaled$se_f,2*z$se_f))

}
# The full working correlation matrix independently defines admissibility.
grid <- expand.grid(a=seq(-1,1,length.out=21),b=seq(-1,1,length.out=21))
positive <- vapply(seq_len(nrow(grid)),function(i) {
  a <- grid$a[i]; b <- grid$b[i]
  R <- matrix(c(1,rxy,rxz,a,rxy,1,ryz,b,rxz,ryz,1,0,a,b,0,1),4)
  min(eigen(R,symmetric=TRUE,only.values=TRUE)$values)>1e-12
},logical(1))
z <- pointwise(grid$a,grid$b,sy/sx,delta_min)
stopifnot(identical(z$valid,positive),all(is.finite(z$t[z$valid])),
          all(is.finite(z$p_delta[z$valid])),all(is.na(z$t[!z$valid])),
          all(is.na(z$p_zero[!z$valid])),all(is.na(z$p_delta[!z$valid])))
# Approach D=0 from both sides: finite valid probabilities, excluded outside.
edge <- sqrt((1-ryz^2)*(1-.13^2))
z <- pointwise(rep(0,2),edge*c(1-1e-12,1+1e-12),sy/sx,delta_min)
stopifnot(identical(z$valid,c(TRUE,FALSE)),is.finite(z$t[1]),
          is.finite(z$p_delta[1]),is.na(z$p_delta[2]))
varying <- pointwise(c(0,.3),c(0,-.3),sy/sx,delta_min)
stopifnot(diff(varying$se_f)!=0,diff(varying$d_delta)!=0,
          all(is.finite(z$se_f[z$valid])),all(is.na(z$se_f[!z$valid])))
cat('PASS: Delta=2.706; A10 direct OLS and T identity; varying SE/cutoff; no extra significance condition; scale invariance; admissibility.\n')
