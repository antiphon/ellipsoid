# Test the uniform ellipse

library(devtools)
load_all(".")

axes <- c(1,5)
Rprof()
x <- rellipsoid(n<-5000, axes)
xf <- rellipsoid_dev(n, axes, method=2)
xfr <- rellipsoid_dev(n, axes, method=1, rej=T)
Rprof("")
mem <- summaryRprof()



diag<-function(x,...){
  plot(x, asp=1, pch=".", ...)
  #hist(x[,1], ...)
  #hist(x[,2], ...)
  pieces <- 100
  theta <- seq(0, 2 * pi, length=pieces)
  xy <- cbind(axes[1] * cos(theta), axes[2] * sin(theta))
  # the surface integrals ~ distance between consecutive points
  integ <- rowSums((xy[-1,] - xy[-pieces,])^2)
  den <- integ/sum(integ)
  # bin
  thet <- atan2(x[,2],x[,1]) # angles of data
  thet[thet<0] <- 2*pi + thet[thet<0]
  z <- cut(thet, theta, labels=FALSE)
  ncount <- tabulate(z, nbins = pieces-1)
  # the number of points should inversely related to the size of the surface
  plot(den, 1/ncount)
  # the sample is better
}





par(mfrow=c(3,2))
diag(x, main="sample")
diag(xf, main="rej")
diag(xfr, main="biased")
# should be a line...
## Sample is better than rejection: Making sample the default.

