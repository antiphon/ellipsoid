#' Project to a plane, 2d outcome
#' 
#' Plane defined using a point and a normal
#' 
#' @param xyz points to be projected, xyz-matrix
#' @param normal normal of the plane
#' 
#' We will simply rotate the points xyz and drop z-coordinate.
#' Rotation order: 1) Rz 2) Ry
#' 
#' This is not unique w.r.t. rotation around the normal, be aware.
#' 
#' @export

project_to_plane <- function(xyz, normal){
  normal <- normal / sqrt(sum(normal^2))
  #'
  a <- atan2(normal[2], normal[1])
  i <- acos(normal[3])
  #' all we need is to straighten the z coordinate:
  Az <- matrix(c(cos(a), -sin(a), 0, sin(a), cos(a),0,0,0,1), 3)
  Ay <- matrix(c(cos(i), 0, sin(i), 0,1,0, -sin(i),0,cos(i)), 3)
  Ay%*%Az%*%normal
  #'
  xy <-t( Ay%*%Az%*%t(xyz) )[,1:2]
  attr(xy,"inc") <- i
  attr(xy,"azi") <- a
  xy
}

#' Three points to a plane
#' 
#' @export
three_points_to_plane <- function(x){
  x <- rbind(x)
  if(nrow(x)!=3) stop("x should be 3 non-collinear points.")
  # normal
  n <- cross(x[2,]-x[1,], x[3,]-x[1,])
  n <- n/sqrt(sum(n^2))
  # location point, take geometric mean
  p <- colMeans(x)
  c(normal=n, p=p)
}


#' Cross-product of two 3D vectors
#' 
#' @export
cross <- function(u, v){
  c(u[2]*v[3]-u[3]*v[2], 
    u[3]*v[1]-u[1]*v[3],
    u[1]*v[2]-u[2]*v[1])
}

