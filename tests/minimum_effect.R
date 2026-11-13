# Run from repository root: Rscript tests/minimum_effect.R
# Load the real evaluator and cutoff without running the Monte Carlo driver.
exprs <- as.list(parse('minimum_effect.R'))
for (name in c('pointwise','d_fixed')) {
  definition <- Filter(function(e) is.call(e) && identical(e[[1]],as.name('<-')) &&
                        identical(e[[2]],as.name(name)),exprs)
  if(length(definition)!=1L) stop('Missing fixed-cutoff definition: ',name)
  eval(definition[[1]])
}
options(warn=2)
stopifnot(abs(d_fixed-2.706)<1e-14)

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
  fit <- lm(y~x+z+c,data=dat)
  tt <- coef(summary(fit))['x','t value']
  z <- pointwise(a,b,d_fixed)
  stopifnot(z$valid,df.residual(fit)==7635,abs(z$t-tt)<1e-10,
            abs(z$p_zero-pnorm(-tt-1.96))<1e-10,
            abs(z$p_fixed-pnorm(-tt-2.706))<1e-10,
            z$p_fixed<=z$p_zero)
  # Changing only the cutoff to 1.96 must reproduce NHST exactly.
  original <- pointwise(a,b,1.96)
  stopifnot(identical(z$t,original$t),identical(z$valid,original$valid),
            identical(original$p_fixed,original$p_zero))
  # At -tt=cutoff the normal argument is zero, hence probability 0.5.
  stopifnot(abs(pointwise(a,b,-tt)$p_fixed-.5)<1e-10)
}
# The full working correlation matrix independently defines admissibility.
grid <- expand.grid(a=seq(-1,1,length.out=21),b=seq(-1,1,length.out=21))
positive <- vapply(seq_len(nrow(grid)),function(i) {
  a <- grid$a[i]; b <- grid$b[i]
  R <- matrix(c(1,rxy,rxz,a,rxy,1,ryz,b,rxz,ryz,1,0,a,b,0,1),4)
  min(eigen(R,symmetric=TRUE,only.values=TRUE)$values)>1e-12
},logical(1))
z <- pointwise(grid$a,grid$b,d_fixed)
stopifnot(identical(z$valid,positive),all(is.finite(z$t[z$valid])),
          all(is.finite(z$p_fixed[z$valid])),all(is.na(z$t[!z$valid])),
          all(is.na(z$p_zero[!z$valid])),all(is.na(z$p_fixed[!z$valid])))
# Approach D=0 from both sides: finite valid probabilities, excluded outside.
edge <- sqrt((1-ryz^2)*(1-.13^2))
z <- pointwise(rep(0,2),edge*c(1-1e-12,1+1e-12),d_fixed)
stopifnot(identical(z$valid,c(TRUE,FALSE)),is.finite(z$t[1]),
          is.finite(z$p_fixed[1]),is.na(z$p_fixed[2]))
cat('PASS: cutoff 2.706; unchanged OLS statistic; cutoff-only substitution; NHST identity; admissibility and boundary stability.\n')
