# from sphere, to reduce dependencies



#' (azi,incl) to 3d coordinates
#' 
#' @export
ai2xyz <- function(aziinc) {
  aziinc <- rbind(aziinc)
  azi <- aziinc[,1]
  inc <- aziinc[,2]
  wx<-sin(inc)*cos(azi)
  wy<-sin(inc)*sin(azi)
  wz<-cos(inc)
  cbind(x=wx,y=wy,z=wz)
}


#' Sample unit sphere uniformly
#' 
#' @param n sample size
#' @param spherical Return azimuth-inclination? (Def: FALSE)
#' @details
#' Default is to return unit vectors.
#' @export

runifsphere <- function(n, spherical=FALSE){
  u <- runif(n)
  v <- runif(n)
  a <- u * 2 * pi
  i <- acos(2*v - 1)
  ai<- cbind(azi=a, inc=i)
  if(spherical) ai
  else ai2xyz(ai)
}


#' Sample unit circle uniformly
#' 
#' @param n sample size
#' @param polar Return polar coordinates? (Def: FALSE)
#' @details
#' Default is to return unit vectors.
#' @export

runifcircle <- function(n, polar=FALSE){
  a <- runif(n,  0, 2 * pi)
  if(polar) cbind(r=1, phi=a)
  else cbind(cos(a),sin(a))
}


#' Product of elementary rotation matrices
#' 
#' @param order Order in which to to the product. Aimed for pre-product
#' Default: 3-2-1 i.e. around z, then around y, then around x. 
#' @export 
rotationMatrix3 <- function(ax=0, ay=0, az=0, order=c(3,2,1)) {
  R<-list()
  R[[1]] <- cbind(c(1,0,0), c(0, cos(ax), sin(ax)), c(0, -sin(ax), cos(ax)) )
  R[[2]] <- cbind(c(cos(ay), 0, -sin(ay)), c(0, 1, 0), c(sin(ay), 0, cos(ay)) )
  R[[3]] <- cbind(c(cos(az), sin(az), 0), c(-sin(az), cos(az), 0), c(0,0,1) )
  R[[order[3]]]%*%R[[order[2]]]%*%R[[order[1]]]
}



#' Rotation matrix to Euler angles
#' 
#' @return
#' (Heading, Attitude, Bank)
#' 
#' @export
rotationMatrix2EulerAngles <- function(R){ 
  quaternion2EulerAngles(rotationMatrix2quaternion(R))
}

#' Quaternion to Euler angles
#' 
#' from euclideanspace.com,
#' 
#' @details
#' 
#' Convention: Euler is (Heading, Attitude, Bank) so that
#' Heading = rotation  around y, Attitude = rotation around z, 
#' Bank = rotation around x, in that order.
#' 
#' @export
quaternion2EulerAngles <- function(q){
  #   if(q[2]*q[3]+q[3]*q[1]==0.5){
  #     phi <- 2*atan2(q[2],q[1])
  #     
  #   }
  phi <- atan2(2*(q[1]*q[3]-q[2]*q[4]), 1-2*(q[3]^2+q[4]^2))
  theta <- asin(2*(q[2]*q[3]+q[1]*q[4]))
  psi <- atan2(2*(q[1]*q[2]-q[3]*q[4]), 1-2*(q[2]^2+q[4]^2))
  c(phi=phi, theta=theta, psi=psi)
}


#' Rotation matrix to quaternion
#' 
#' @export
rotationMatrix2quaternion <- function(R){
  q4 <- 0.5 * sqrt(1 + R[1,1] + R[2,2] + R[3,3])
  if(is.na(q4)) q4<-0 # improper rotation
  if(q4 < 1e-5){ # try other way
    q1 <- 0.5 * sqrt(1+R[1,1]-R[2,2]-R[3,3])
    if(q1 < 1e-5){ # try other way
      q2 <- 0.5 * sqrt(1-R[1,1]+R[2,2]-R[3,3])
      if(q2 < 1e-5){
        q3 <- 0.5 * sqrt(1-R[1,1]-R[2,2]+R[3,3])
        if(q3 < 1e-5) stop("Can't convert matrix to quaternion.")
        m <- 1/(4*q3) # q3 ok
        q1 <- m * (R[3,1]+R[1,3])
        q2 <- m * (R[3,2]+R[2,3])
        q4 <- m * (R[2,1]-R[1,2])
      }else{ # q2 ok
        m <- 1/(4*q2)
        q1 <- m * (R[2,1]+R[1,2])
        q3 <- m * (R[3,2]+R[2,3])
        q4 <- m * (R[1,3]-R[3,1])
      }
    }else{ # q1 ok 
      m <- 1/(4*q1)
      q2 <- m * (R[1,2]+R[2,1])
      q3 <- m * (R[1,3]+R[3,1])
      q4 <- m * (R[3,2]-R[2,3])
    }
  }else{ # q4 ok
    m <- 1/(4*q4)
    q1 <- m * (R[3,2]-R[2,3])
    q2 <- m * (R[1,3]-R[3,1])
    q3 <- m * (R[2,1]-R[1,2])
  }
  c(q1,q2,q3,q4)
}


#' values to colors
#' 
#' @param v the values
#' @param n number of colors
#' @param zlim limits
#' @param color function, e.g. heat.colors, gray.colors
#' 
#' @export

values2colors <- function(v, n=100, zlim, col=heat.colors, na.col="gray50", ...){
  zlim0 <- range(v, na.rm = TRUE)
  if(missing(zlim)) zlim <- zlim0
  pr <- v
  na <- which(is.na(pr))
  pr[na] <- mean(pr, na.rm=TRUE)
  pr <- pmax(zlim[1], pmin(zlim[2], pr))
  pp <- (pr-zlim[1])/diff(zlim)
  e <- floor(pp*(n-1)) + 1
  hc <-col(n)
  out <- hc[e]
  out[na] <- na.col
  out
}

