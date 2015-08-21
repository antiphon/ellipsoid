############################################################
#' Plot an ellipsoid
#' 
#' @exportMethod plot
#' @export
plot.ellipsoid <- function(x, add=TRUE, res=201, scale=1, ...){
  if(x$dim==2){
    a <- c(seq(0, 2*pi, length=res))
    y <- cbind(cos(a),sin(a))
    z <- y * predict(x, y) * scale
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
#'@exportMethod predict
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
#' @exportMethod print
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
#' @export
ellipse_form <- function(beta, d){
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
  # check definitenes
  e <- eigen(Ahat)
  if(any(e$values < 0)){
    i <- which(e$values > 0)
    Av <- diag(0, ncol(A))
    for(j in i) Av <- Av + e$val[j] * e$vec[,j]%*%t(e$vec[,j])
    Ahat <- Av
  }
  list(c=chat, A=Ahat)
}

