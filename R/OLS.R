#############################################################
#' OLS fitting of an ellipsoid to a d-dimensional point cloud
#'
#' @param x point coordinate matrix
#' @param origin Fix center of the ellipsoid to be the origin?
#' 
#' @export
ellipsoid_OLS  <- function(x, origin=FALSE) {
  n <- nrow(x)
  d <- ncol(x)
  if(n < 3+d) stop("Trying to fit an ellipsoid to too little amount of points.")
  D  <- matrix(2, ncol=d, nrow=d)
  diag(D) <- 1
  D <- D[upper.tri(D,T)]
  ## Y:
  Y<-apply(x, 1, function(x){
    X  <- x%*%t(x)
    vX <- X[upper.tri(X, T)]
    y <- c(D*vX, if(origin) NULL else c(x), 1)
    y
  })
  Y <- t(Y)
  #
  nd <- (d*(d+1)/2)
  nb <- nd + (d * !origin) + 1 
  H <- diag(1, nb)
  Hi <- solve(H)
  # svd
  USV <- svd(Y%*%Hi)
  # the estimate of parameters:
  beta <- c( Hi%*%USV$v[,nb] )
  if(origin){ # add center parameters to comply with other functions
    beta <- c(beta[1:nd], rep(0, d), beta[(nd+1):nb])
  }
  
  # convert to center and matrix -formulation
  res <- ellipsoid_from_beta(beta, d)
  res$ndata <- n
  ####### add the error variance
  # the error variance
  r_data <- sqrt(rowSums(x^2))
  u_data <- x/r_data
  r_pred <- predict(res, u_data)
  resi <- r_data - r_pred
  if(d==2){
    s2 <- sum(resi^2)/(n-1) # this is approximation with 0 curvature TODO better
  }
  else{
    s2 <- sum(resi^2)/(n-1) # this is approximation to 0 curvature TODO BETTER
  }
  # the variance covariance
  warn <- FALSE; m <- 1e-9 # hedging
  ss <- try(S0 <- solve(t(Y)%*%Y + diag(m, ncol(Y))), TRUE)
  while("try-error"%in% is(ss)){
    ss <- try(S0 <- solve(t(Y)%*%Y + diag(m <- m * 2, ncol(Y))), TRUE)
    warn <- TRUE
  }
  if(warn) warning(paste("Hedging needed for solving covariance matrix. (diag ", m, ")"))
  S <- s2 * S0
  #
  if(origin){ # add dummy for center
    S0 <- diag(1e-9, nb+d)
    S0[1:(nb-1), 1:(nb-1)] <- S[-nb,-nb]
    aa <- c(1:(nb-1), nb+d)
    S0[aa, nb+d] <- S0[nb+d, aa] <- S[nb,]
    S0 <- 0.5 * (S0 + t(S0))
    S <- S0
  }
  # make sure symmetric varcov
  
  if(max(S-t(S))>0.01) warning("varcov not good, big asymmetry.")
  S <- 0.5 * (S + t(S))
  #
  res$ols_fit <- list(varcov=S, beta_est=beta, s2=s2)
  res$origin <- origin
  #
  #' done
  res
}

#########################################################################
#' Solve rotation and semi-axes from general transform
#' @param A the trasformation matrix in ellipsoid equation
#' @export

ellipse_solve_rota <- function(A){
  ev <- eigen(A)
  # make sure working with a definite:
  if(any(ev$value<0)){
    i <- which(ev$value > 0)
    S <- diag(0, ncol(A))
    for(j in i) S <- S + ev$value[j] * ev$vec[,i]%*%t(ev$vec[,i])
    ev <- eigen(S)
  }
  axes_len <- 1/sqrt(ev$value) # the semi-axes lengths
  R <- ev$vector # this holds the rotation...
  # check proper
  #   if(det(R)<0){
  #     # mirror
  #     R <- Euler2rotationMatrix(rotationMatrix2Euler(R))
  #   }
  list(axes=axes_len, R=R)
}


##########################################################################
#' Sample from the OLS estimate of the beta parameters
#'
#' @import mvtnorm 
#' @export
sample_ellipse_beta <- function(x, nsim=100, tol=0, maxiter=5000){
  d <- x$dim
  nb  <- (d*(d+1)/2) + d + 1
  
  b <- rmvnorm(nsim, x$ols_fit$beta_est, x$ols_fit$varcov)
  if(maxiter==0){ tol<-0; b <- b/sqrt(rowSums(b^2))}
  if(tol>0){
    dev <-  abs(sqrt(rowSums(b^2))-1)
    ok <- dev < tol
    b <- b[ok,]
    it <- 0
    failed <- FALSE
    while(length(b)/nb < nsim){
      b <- rbind(b, rmvnorm( 2*(nsim-length(b)/nb), x$ols_fit$beta_est, x$ols_fit$varcov))
      ok <- abs(sqrt(rowSums(b^2))-1) < tol
      b <- b[ok,]
      it<-it+1
      if(it>maxiter) {it<-0; tol <- tol * 10; failed<-TRUE}
    }
    if(failed) warning(paste("beta sampling tolerance was increased to", tol))
    b<-b[1:nsim,]
  }
  b
}


#####################################################################################
#' convert beta vector to ellipsoid object
#' 
#' @param beta beta vector, the coefficients in quadratic form
#' @param d dimension
#' @import sphere
#' @export

ellipsoid_from_beta <- function(beta, d){
  elform <- ellipse_form(beta, d)
  chat <- elform$c
  Ahat <- elform$A
  # solve rotation and axes:
  rota <- ellipse_solve_rota(Ahat)
  R <- rota$R
  axes_len <- rota$axes
  M <- R%*%diag(axes_len)
  angles <- NULL
  if(d==2) {
    f <- R %*% c(1,0)
    angles <- atan2(f[2],f[1])
  }else if(d==3){
    angles <- sphere::rotationMatrix2EulerAngles(R)
  }
  
  # check if we got a valid fit
  valid <- all(!is.infinite(c(axes_len, angles)) & !is.na(c(axes_len, angles)))
  # compile
  
  res <- list(center=chat, A=Ahat, 
              semi_axes=axes_len, 
              rot=R, 
              M=M, 
              rot_angle=angles, 
              valid=valid, dim=d)
  class(res) <- "ellipsoid"
  res
}

