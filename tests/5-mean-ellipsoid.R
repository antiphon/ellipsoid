#' test the mean ellipsoid
library(devtools)
load_all(".")

# some ellipsoids
els <- list()
sds <- seq(0.05, 0.25, length=5)
g <- runif(5, 0.5, 1.5)
for(i in 1:5) els[[i]] <- rellipsoid(rpois(1,40), axes=c(g[i],1/g[i]), noise.sd = sds[i], center=c(0,1))

N <- sum(sapply(els, nrow))
# fit ellipses
fits <- lapply(els, ellipsoid_OLS, origin = F)

m <- mean_ellipsoids(fits, keep_data = TRUE, add_noise = TRUE, nsim = N, origin = F)
m2 <- mean_ellipsoids(fits, keep_data = TRUE, add_noise = TRUE, nsim = N*5, origin = F)
m3 <- mean_ellipse(fits)
#

par(mfrow=c(1,2), mar=c(0,0,0,0))
zc <- 2
plot(els[[1]], col=1, asp=1, xlim=c(-1,1)*zc, ylim=c(-1,1)*zc)
for(i in 1:5) {
  points(els[[i]], col=i, pch=i)
  plot(fits[[i]], col=i)
}
legend("topleft", legend=1:5, fill=1:5)

plot(els[[1]], asp=1,xlim = c(-1,1)*zc, ylim=c(-1,1)*zc)
plot(m, col=1, lwd=5)
points(m$data, col=1, cex=.2)
plot(m2, col=2, lwd=5)
points(m2$data, col=2, pch=2, cex=.2)
plot(m3, col=3, lwd=5)

#points(m$data, pch=19, cex=0.4)

cat("points :", sapply(els, nrow))
cat("\nestd sd:", round(sqrt(s2<-sapply(lapply(fits, getElement, "ols_fit"), getElement, "s2")),2)  )
cat("\ntrue sd:", round(sds,2))
cat("\n order :", order(s2))
cat("\n m sims:", (m$nsims))
