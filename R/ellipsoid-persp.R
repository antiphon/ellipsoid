#' Perspective plot of an ellipsoid in 3d
#' 
#' Try to make a pretty perspective plot of an 3D ellipsoid.
#' 
#' @param x ellipsoid-object
#' @param add Add to a scene? 
#' @param theta Camera rotation, see persp.default
#' @param phi Camera rotation, see persp.default
#' @param expand zoom 
#' @param triptych plot from all three main axes?
#' @param pmat optional projection matrix (e.g. from 2d persp)
#'
#' @import graphics
#' @export

persp.ellipsoid <- function(x, add=FALSE, theta=25, phi=30, 
                            expand=.9, triptych=FALSE, 
                            pmat, ...){
  if(x$dim != 3) stop("perspective plot only for 3D ellipsoids.")
  if(triptych){
   a <- c(90, 180, 0)
   b <- c(0, 0,90)
   for(i in 1:3) pmat <- persp.ellipsoid(x, add, theta=a[i], 
                                         phi=b[i], expand, triptych=FALSE, pmat, ...)
  }
  else{
    L <- c(-1,1)*max(x$semi_axes) * 1
    if(!add){
      pmat <- persp(L, L, matrix(NA, 2, 2), 
                  xlim=L, ylim=L, zlim=L,
                  theta=theta, phi=phi, expand=expand , 
                  xlab="X", ylab="Y", zlab="Z", 
                  box=TRUE, ticktype="simple", shade=FALSE)
    }else 
      if(missing(pmat)) 
        stop("Need to provide 'pmat'-projection matrix for adding to persp-plot.")
    #'
    add_ellipsoid2persp(x=x, pmat=pmat, ...)
    #
  }
  invisible(pmat)
}

#' add ellipse to persp plot
#' 
#'
add_ellipsoid2persp <- function(x, N=2, colmap=values2colors, pmat, ...){
  s <- ellipsoid_shape(N, x$semi_axes, x$rot)
  w <- s$vb[4,]
  xc <- s$vb[1,]/w
  yc <- s$vb[2,]/w
  zc <- s$vb[3,]/w
  coords <- cbind(xc,yc,zc)
  Lm <- range(zc)
  zm <- apply(s$it, 2, function(n) mean(zc[n]) )
  cols <- colmap(zm)
  # ah the problem is plotting of behind after front:
  # order according to camera position:
  camera <- c(unlist(trans3d(0,0,0, pmat)),100)
  # distance surface
  mini <- apply(s$it, 2, function(n) 
  {d<-t(t(coords[n,]) - camera); n[which.min(diag(d%*%t(d)))]}  )
  o <- order(coords[mini,3])
  for (j in o) {
    idx <- s$it[,j]
    co <- coords[idx,]
    col <- cols[j]
    polygon(trans3d(co[,1], co[,2], co[,3], pmat), col=col, ...)
  }
}
