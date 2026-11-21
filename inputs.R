# Shared working inputs. Education is computed from Frank (2000), Table 1.
# Retention remains the explicitly assumption-dependent rounded reference.
education_inputs <- local({
  rxy <- 0.427
  rXZ <- c(-0.437, -0.289)
  rYZ <- c(-0.317, -0.350)
  RZZ <- matrix(c(1, 0.263, 0.263, 1), 2L)
  cross_term <- drop(crossprod(rXZ, solve(RZZ, rYZ)))
  rxz2 <- drop(crossprod(rXZ, solve(RZZ, rXZ)))
  ryz2 <- drop(crossprod(rYZ, solve(RZZ, rYZ)))
  R <- rbind(c(1, rxy, rXZ), c(rxy, 1, rYZ), cbind(rXZ, rYZ, RZZ))
  list(partial_r=(rxy-cross_term)/sqrt((1-rxz2)*(1-ryz2)),
       rxz2=rxz2, ryz2=ryz2, nu=10562, direction=1,
       cross_term=cross_term, correlation_matrix=R)
})
retention_inputs <- list(partial_r=-0.13, rxz2=0.2704, ryz2=0.2025,
                         nu=7635, direction=-1, sx=0.240548, sy=13.53)
