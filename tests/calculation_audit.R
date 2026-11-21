# Run from the repository root: Rscript tests/calculation_audit.R
# Definitions are loaded without running either Monte Carlo driver.
options(warn=2)
source('inputs.R')
close <- function(x,y,tol=1e-12) {
  stopifnot(isTRUE(all.equal(x,y,tolerance=tol,check.attributes=FALSE)))
}
exprs <- as.list(parse('robustness.R'))
for (name in c('evaluate_pairs','make_surface')) {
  definition <- Filter(function(e) is.call(e) && identical(e[[1]],as.name('<-')) &&
                         identical(e[[2]],as.name(name)),exprs)
  stopifnot(length(definition)==1L)
  eval(definition[[1]])
}
threshold <- .8
# Independently specified Table 1 correlation matrix (x,y,z1,z2).
R_edu <- matrix(c(1,.427,-.437,-.289,
                  .427,1,-.317,-.350,
                  -.437,-.317,1,.263,
                  -.289,-.350,.263,1),4,byrow=TRUE)
stopifnot(min(eigen(R_edu,symmetric=TRUE,only.values=TRUE)$values)>0)
close(education_inputs$correlation_matrix,R_edu)
close(education_inputs$cross_term,.18838965505016483)
close(education_inputs$rxz2,.22352057677494624)
close(education_inputs$ryz2,.17686271729239786)
close(education_inputs$partial_r,.2984613894645406)
stopifnot(education_inputs$nu==10562,education_inputs$direction==1)
# Precision-matrix partial correlation, independent of the production formula.
precision <- solve(R_edu)
r_precision <- -precision[1,2]/sqrt(precision[1,1]*precision[2,2])
close(education_inputs$partial_r,r_precision)
# FWL via explicit projection on two standardized controls.
n <- 10567L
k <- seq_len(n)
basis <- qr.Q(qr(cbind(1,sin(k),cos(k),sin(2*k),cos(2*k))))[,2:5]*sqrt(n-1)
data_edu <- basis %*% chol(R_edu)
controls <- cbind(1,data_edu[,3:4])
x_resid <- qr.resid(qr(controls),data_edu[,1])
y_resid <- qr.resid(qr(controls),data_edu[,2])
close(cor(x_resid,y_resid),r_precision)
T0 <- r_precision*sqrt(10562/(1-r_precision^2))
close(T0,32.13815,tol=2e-7)
close(retention_inputs$partial_r,-.13)
close(retention_inputs$rxz2,.2704); close(retention_inputs$ryz2,.2025)
close(retention_inputs$sx,.240548); close(retention_inputs$sy,13.53)
stopifnot(retention_inputs$nu==7635,retention_inputs$direction== -1)
rxz <- .52; ryz <- -.45
rxy <- -.13*sqrt((1-rxz^2)*(1-ryz^2))+rxz*ryz
R_ret <- matrix(c(1,rxy,rxz,rxy,1,ryz,rxz,ryz,1),3)

for (example in c('Retention','Education')) {
  s <- if(example=='Education') education_inputs else retention_inputs
  R <- if(example=='Education') R_edu else R_ret
  for (bound in c(1,.5)) {
    axes <- seq(-bound,bound,length.out=21)
    grid <- expand.grid(a=axes,b=axes)
    # Schur complement of the joint control block gives residual covariance V.
    # The full correlation matrix includes c independent of the original Z.
    point <- lapply(seq_len(nrow(grid)),function(i) {
      a <- grid$a[i]; b <- grid$b[i]
      full <- rbind(cbind(R,c(a,b,rep(0,nrow(R)-2))),
                    c(a,b,rep(0,nrow(R)-2),1))
      controls <- 3:nrow(full)
      V <- full[1:2,1:2] - full[1:2,controls,drop=FALSE] %*%
        solve(full[controls,controls],full[controls,1:2,drop=FALSE])
      valid <- min(eigen(full,symmetric=TRUE,only.values=TRUE)$values)>1e-12
      if(valid) {
        stopifnot(min(eigen(V,symmetric=TRUE,only.values=TRUE)$values)>0)
        tt <- V[1,2]*sqrt(s$nu/det(V))
        p <- pnorm(s$direction*tt-1.96)
        stopifnot(is.finite(tt),is.finite(p),p>=0,p<=1)
      } else p <- NA_real_
      c(valid=valid,p=p,D=det(V),vx=V[1,1],vy=V[2,2])
    })
    oracle <- do.call(rbind,point)
    valid <- as.logical(oracle[,'valid'])
    result <- evaluate_pairs(axes,axes,s$partial_r,s$rxz2,s$ryz2,s$nu,s$direction)
    stopifnot(result['n_nonfinite']==0,result['n_admissible']==sum(valid))
    close(result['n_above'],sum(oracle[valid,'p']>.8))
    close(result['p_hat'],mean(oracle[valid,'p']>.8))
    close(result['min_D'],min(oracle[valid,'D']))
    close(result['min_vx'],min(oracle[valid,'vx']))
    close(result['min_vy'],min(oracle[valid,'vy']))
    surface <- make_surface(-bound,bound,s$partial_r,s$rxz2,s$ryz2,
                            s$nu,s$direction,n_grid=21L)
    stopifnot(identical(as.vector(surface$admissible),valid),
              all(is.na(surface$z[!surface$admissible])),
              all(is.finite(surface$z[surface$admissible])))
    close(as.vector(surface$z)[valid],oracle[valid,'p'])
  }
  # Direct OLS with exact sample covariance tests the t statistic and its sign.
  n <- s$nu+nrow(R)+1L
  k <- seq_len(n)
  basis <- qr.Q(qr(cbind(1,sin(k),cos(k),sin(2*k),cos(2*k),sin(3*k))))
  basis <- basis[,2:(nrow(R)+2L)]*sqrt(n-1)
  for (ab in list(c(0,0),c(.3,-.3),c(-.4,-.3),c(.5,.46),c(.4,-.23))) {
    full <- rbind(cbind(R,c(ab,rep(0,nrow(R)-2))),
                  c(ab,rep(0,nrow(R)-2),1))
    stopifnot(min(eigen(full,symmetric=TRUE,only.values=TRUE)$values)>0)
    dat <- as.data.frame(basis %*% chol(full))
    names(dat) <- c('x','y',paste0('z',seq_len(nrow(R)-2)),'c')
    fit <- lm(reformulate(setdiff(names(dat),'y'),response='y'),dat)
    tt <- coef(summary(fit))['x','t value']
    stopifnot(df.residual(fit)==s$nu,is.finite(tt))
    result <- evaluate_pairs(ab[1],ab[2],s$partial_r,s$rxz2,s$ryz2,s$nu,s$direction)
    close(result['mean_pointwise_pmc'],pnorm(s$direction*tt-1.96))
    close(result['p_hat'],as.numeric(s$direction*tt>qnorm(.8)+1.96))
  }
}
cat('PASS: Education Table 1, precision matrix and FWL; retained Retention inputs; all four scenario grids, full-matrix/V positive definiteness, finite probabilities and direct OLS.\n')
