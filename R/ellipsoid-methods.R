############################################################
#' Plot an ellipsoid
#' 
#' @export
plot.ellipsoid <- function(x, add=TRUE, res=201, scale=1, ...){
  if(x$dim==2){
    a <- c(seq(0, 2*pi, length=res))
    y <- cbind(cos(a),sin(a))
    z <- y * predict(x, y) * scale
    if(!add) plot(NA, xlim=range(z), ylim=range(z), asp=1, xlab="", ylab="")
    lines(z, ...)
  }
  else{
      y <- ellipsoid_shape(axes=x$semi_axes, R=x$rot)
      rgl::shade3d(y, ...)
  }
}
####################################################################
#' Predict i.e. give the length of a direction to be on the ellipsoid
#'
#' Returns the distance from ellipsoid center to the ellipsoid surface in
#' the given directions.
#'
#'@export
predict.ellipsoid <- function(x, u, ...){
  if(missing(u)) stop("direction(s) u needed")
  d <- 1 / diag(u %*% x$A %*% t(u) )
  r <- sqrt(d) 
  r
}



####################################################################
#' Ellipsoid shape for 3d plotting
#'
#' Refine an icosahedron, then transform
#' 
#' @param R Rotation matrix.
#'
#' @import rgl
#' @export

ellipsoid_shape <- function(N=2, axes=c(1,1,1), R=NULL){
  ico <- rgl::icosahedron3d()
  for(i in 1:N) ico <- rgl::subdivision3d(ico)
  D <- diag(axes)
  xy <- t(D%*%ico$vb[1:3,])
  exy <- if(is.null(R)) xy else t(R%*%t(xy))
  ico$vb[1:3,] <- t(exy)
  ico$vb <- t( t(ico$vb)/apply(ico$vb, 2, function(v) sqrt(sum(v^2))))
  ico
}


############################################################
#' Print ellipsoid
#' 
#' @export
print.ellipsoid <- function(x, ...){
  type <- ifelse(x$dim==2, "2D ellipse", "3D ellipsoid")
  if(!is.null(x$ave))
    cat(paste0("Average ", type, ", computed from ", x$nellipses, " ", type, "s.\n"))
  else cat(type, "fitted to", x$n, "points.\n")
}

####################################################################
#' Ellipse center and matrix from general parameter form
#' 
#' @param beta OLS estimates
#' @param d dimension
#' @param check Check for definiteness?
#' 
#' @export
ellipse_form <- function(beta, d, check=FALSE){
  nd <- (d*(d+1)/2)
  nb <- nd + d + 1 
  # ellipse parameters
  A <- diag(0, d)
  A[upper.tri(A,T)] <- beta[1:nd]  
  A[lower.tri(A)] <- A[upper.tri(A)]
  b <- beta[(nd+1):(nd+d)]
  dhat <- beta[nb]
  m <- 0
  v <- try(Sa <- solve(A + diag(m, ncol(A))))
  while("try-error"%in% is(v))  v <- try(Sa <- solve(A + diag(m<-m+5e-8, ncol(A))))
  chat <- -0.5 * Sa%*%b
  Ahat <- A / (t(chat) %*% A %*% chat - dhat)[1]
  # check definiteness
  if(check){
    e <- eigen(Ahat)
    if(any(e$values < 0)){
      i <- which(e$values > 0)
      Av <- diag(0, ncol(Ahat))
      for(j in i) Av <- Av + e$val[j] * e$vec[,j]%*%t(e$vec[,j])
      Ahat <- Av
    }
  }
  list(c=chat, A=Ahat)
}

#########################################################################
#' Solve rotation and semi-axes from general transform
#' 
#' @param A the trasformation matrix in an ellipsoid equation
#' @param eps Inflate diagonal by this factor, to avoid numerical problems.
#' 
#' @export

ellipse_solve_rota <- function(A, eps = 0){
  ev <- eigen(A + diag(diag(A)*eps))
  # make sure working with a definite:
  if(any(ev$value<0)){
    i <- which(ev$value > 0)
    S <- diag(0, ncol(A))
    for(j in i) S <- S + ev$value[j] * ev$vec[,i]%*%t(ev$vec[,i])
    ev <- eigen(S)
  }
  # in case negative eigenvalues persist, they are
  # extremely close to -0 so we might as well take abs here.
  axes_len <- 1/sqrt(abs(ev$value)) # the semi-axes lengths
  R <- ev$vector # rotation
  list(axes=axes_len, R=R)
}


##########################################################################
#' Sample from the OLS estimate of the beta parameters
#' 
#' Simulate the beta parameters, approximately normal conditional on ||beta||^2=1
#' 
#' @param x ellipsoid object
#' @param nsim number of simulations
#' @param tol tolerance when comparing to vector length 1
#' @param maxiter maximum number of iterations to try to get enough samples within tolerance
#' @import mvtnorm 
#' @export
sample_ellipse_beta <- function(x, nsim=100, tol=0, maxiter=500){
  d <- x$dim
  nb  <- (d*(d+1)/2) + d + 1
  
  b <- rmvnorm(nsim, x$ols_fit$beta_est, x$ols_fit$varcov)
  if(maxiter==0){ 
    tol <- 0
    b <- b/sqrt(rowSums(b^2))
  }
  if(tol>0){
    dev <-  abs(sqrt(rowSums(b^2))-1)
    ok <- dev < tol
    b <- b[ok,]
    it <- 0
    failed <- FALSE
    while(length(b)/nb < nsim){
      b <- rbind(b, rmvnorm( 2*(nsim-length(b)/nb), 
                             x$ols_fit$beta_est, 
                             x$ols_fit$varcov))
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
#' @param ... passed on to 'ellipse_solve_rota'
#' 
#' @import sphere
#' @export

ellipsoid_from_beta <- function(beta, d, ...){
  elform <- ellipse_form(beta, d)
  chat <- elform$c
  Ahat <- elform$A
  # solve rotation and axes:
  rota <- ellipse_solve_rota(Ahat, ...)
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

