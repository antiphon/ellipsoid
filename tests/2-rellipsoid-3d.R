# Test the uniform ellipse

library(devtools)
load_all(".")
library(rgl)
library(sphere)
library(sp)

axes <- c(1/4,1,4)

x <- rellipsoid_dev(n<-5000, axes)
xf <- rellipsoid(n, axes)


rescale <- function(b) ((b-min(b))/(max(b)-min(b)))

diag<-function(x,...){
  gr <- sphere.tri(3)
  # triangles
  a <- triangulation_areas(gr)
  #
  triangle <- {v<-gr$vb[,gr$it[,1]]; t(v[1:3,])/v[4,]}
  plot3d(gr, aspect=F, col="gold", lit=F)
  points3d(triangle)
  u <- runifsphere(1000)
  # transform
  el <- ellipsoid_shape(3, axes)
  b <- triangulation_areas(el)
  plot3d(el, col=gray(rescale(b)), aspect=F)  
  # shift
  shade3d(translate3d(gr, x=4,0,0), col=gray(rescale(a)))
  plot(a,b)
}
el <- ellipsoid_shape(3, axes)
plot3d(x*.85, aspect=F, size=s<-1)
shade3d(el, lit=F, col="gold")
points3d(t(t(xf*.85)+c(2,0,0)), size=s, col=4)
shade3d(translate3d(el, 2,0,0), lit=F, col="gold")
#plot3d(tri, col=gray(probs[c(tri$it)]))


# should be a line...


