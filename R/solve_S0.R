#' Solve covariance matrix from unstable X
#' 
solve_S0 <- function(Y){
  warn <- FALSE
  m <- 1e-9 # hedging
  ss <- try(S0 <- solve(t(Y)%*%Y + diag(m, ncol(Y))), TRUE)
  while("try-error"%in% is(ss)){
    ss <- try(S0 <- solve(t(Y)%*%Y + diag(m <- m * 2, ncol(Y))), TRUE)
    warn <- TRUE
  }
  if(warn) warning(paste("Hedging needed for solving covariance matrix. (diag ", m, ")"))
  S0
}