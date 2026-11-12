# Run from the repository root: Rscript tests/minimum_effect.R
# Load only the real evaluator; do not run the Monte Carlo driver here.
exprs <- parse('minimum_effect.R')
pointwise_expr <- Filter(function(e) is.call(e) && identical(e[[1]],as.name('<-')) &&
                         identical(e[[2]],as.name('pointwise')),as.list(exprs))
stopifnot(length(pointwise_expr)==1L)
eval(pointwise_expr[[1]])
options(warn=2)

# Independent A10 oracle: exact sample covariance and a direct OLS fit.
n <- 7639L
k <- seq_len(n)
basis <- qr.Q(qr(cbind(1,sin(k),cos(k),sin(2*k),cos(2*k))))[,2:5]*sqrt(n-1)
rxz <- .52; ryz <- -.45
rxy <- -.13*sqrt((1-rxz^2)*(1-ryz^2))+rxz*ryz
sx <- .240548; sy <- 13.374335
for (ab in list(c(0,0),c(.3,-.2),c(-.4,-.3))) {
  a <- ab[1]; b <- ab[2]
  R <- matrix(c(1,rxy,rxz,a,rxy,1,ryz,b,rxz,ryz,1,0,a,b,0,1),4)
  dat <- as.data.frame(basis %*% chol(R))
  names(dat) <- c('x','y','z','c')
  dat$x <- dat$x*sx; dat$y <- dat$y*sy
  fit <- lm(y~x+z+c,data=dat)
  beta <- unname(coef(fit)['x'])
  se <- coef(summary(fit))['x','Std. Error']
  z <- pointwise(a,b,sy/sx,.20*(9.01/.67))
  stopifnot(z$valid,df.residual(fit)==7635,
            abs(z$beta-beta)<1e-10,abs(z$se-se)<1e-10,
            abs(z$t-beta/se)<1e-10)
  # Regression: at the minimum-effect boundary probability is one half.
  at_boundary <- pointwise(a,b,sy/sx,-beta)
  if (abs(at_boundary$p_delta-.5)>1e-10)
    stop('A11 regression: probability at beta=-Delta must be 0.5, without a significance cutoff')
  # Probability below -Delta, independently integrated from the fitted normal.
  expected <- integrate(function(x) dnorm(x,beta,se),-Inf,-.20*(9.01/.67),
                        rel.tol=1e-10)$value
  stopifnot(abs(z$p_delta-expected)<1e-9,
            abs(z$p_zero-pnorm(-beta/se-1.96))<1e-10)
}
# Strict region: invalid residual variances / determinant stay excluded as NA.
z <- pointwise(c(1,0,.8),c(0,1,-.8),sy/sx,.20*(9.01/.67))
stopifnot(!any(z$valid),all(is.na(z$beta)),all(is.na(z$se)),
          all(is.na(z$p_zero)),all(is.na(z$p_delta)))
# Approach D=0 from both sides with a=0: finite probabilities inside only.
edge <- sqrt((1-ryz^2)*(1-.13^2))
z <- pointwise(rep(0,2),edge*c(1-1e-12,1+1e-12),sy/sx,.20*(9.01/.67))
stopifnot(identical(z$valid,c(TRUE,FALSE)),is.finite(z$se[1]),z$se[1]>0,
          is.finite(z$p_delta[1]),is.na(z$p_delta[2]))
cat('PASS: A10 direct OLS; A11 boundary and normal integral; NHST; strict region and boundary stability.\n')
