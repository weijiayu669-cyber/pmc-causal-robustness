# Assumption-dependent reconstruction from Hong & Raudenbush (2005),
# Tables 8 (z; N=7639) and 9 (y; N=7608), p. 216.
# Reference joins Table 8 counts to Table 9 summaries, as in original mcdat.
# outcome_matched uses Table 9 counts for BOTH margins: z summaries in this
# smaller sample are unavailable, so using Table 8 z summaries is an assumption.
# Neither scenario recovers raw records or the unidentified within-cell covariance.
retention_cells <- function(sample=c("reference", "outcome_matched")) {
  sample <- match.arg(sample)
  n <- c(9,12,14,12,24,23,47,48,46,48,47,49,45,39,8,3054,1661,983,323,441,154,213,143,85,49,35,16,8,3)
  meanz <- c(-5.42, -4.25, -3.42, -2.98, -2.57, -2.15, -1.75, -1.24, -0.83, -0.46, 
             -0.01, 0.53, 1.16, 1.99, 3.7, -5.38, -4.27, -3.5, -2.97, -2.55, 
             -2.14, -1.77, -1.27, -0.86, -0.47, -0.04, 0.56, 1.14, 1.8)
  sdz <- c(0.4, 0.24, 0.2, 0.1, 0.17, 0.07, 0.15, 0.14, 0.11, 0.12, 0.11, 0.22, 
           0.2, 0.33, 1.09, 0.39, 0.25, 0.2, 0.1, 0.16, 0.07, 0.15, 0.14, 0.1, 
           0.11, 0.15, 0.17, 0.17, 0.17)
  meany <- c(47.41,47.25,42.66, 45.24,44.07,36.63,37.41,39.79,37.78,35,34.3,32.06,31.09,33.87,32.84,
             59.22,55.62,52.49,52.45,50.05,48.77,45.89,46,43.28,42.18,40.45,33.03,35.27,40.47)
  sdy <- c(11.73,15.50,15.92,17.32,14.43,13.88,10.92,15.37,12.27,9.77,9.05,9.31,11.96,9.15,
           7.71,11.64,11.30,11.95,12.53,12.15,12.67,12.40,13.37,11.52,13.34,13.08,14.33,8.25,8.60)
  if (sample == "outcome_matched") {
    n <- c(9,12,14,12,24,23,47,48,45,48,47,49,45,38,8,
           3044,1654,977,321,440,153,211,143,85,49,35,16,8,3)
  }
  data.frame(n=n, x=c(rep(1,15),rep(0,14)), strata=seq_len(29L),
             meanz=meanz, sdz=sdz, meany=meany, sdy=sdy)
}

mcdat <- function(rho_within=0, sample="reference", return_data=FALSE) {
  if (length(rho_within)!=1L || !is.finite(rho_within) || abs(rho_within)>1) {
    stop("rho_within must be a finite scalar between -1 and 1")
  }
  cells <- retention_cells(sample)
  size <- cells$n
  # Preserve the original sequence: ALL cell-specific z draws, THEN all y draws.
  # At rho=0 use the unmodified original mean/sd calls for bitwise equality.
  vrnorm <- Vectorize(rnorm,c("n","mean","sd"))
  z <- unlist(vrnorm(size,cells$meanz,cells$sdz))
  y <- unlist(vrnorm(size,cells$meany,cells$sdy))
  k <- rep(cells$strata,times=size)
  if (rho_within != 0) {
    ze <- (z-cells$meanz[k])/cells$sdz[k]
    ye <- (y-cells$meany[k])/cells$sdy[k]
    y <- cells$meany[k]+cells$sdy[k]*(rho_within*ze+sqrt(1-rho_within^2)*ye)
  }
  dat <- data.frame(x=rep(cells$x,times=size), strata=k, z=z, y=y)
  if (return_data) return(dat)
  rxz <- cor(dat$x,dat$z)
  ryz <- cor(dat$y,dat$z)
  rxy <- cor(dat$x,dat$y)
  rxyonz <- (rxy-rxz*ryz)/sqrt((1-rxz^2)*(1-ryz^2))
  c(rxz,ryz,rxy,rxyonz)
}
