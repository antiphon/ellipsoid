#' test the mean ellipsoid
library(devtools)
load_all(".")

# some ellipsoids
els <- list()
sds <- seq(0.05, 0.5, length=5)
g <- runif(5, 0.5, 1.5)
for(i in 1:5) els[[i]] <- rellipsoid(rpois(1,40), axes=c(g[i],1/g[i]), noise.sd = sds[i] )

N <- sum(sapply(els, nrow))
# fit ellipses
fits <- lapply(els, ellipsoid_OLS)

m <- mean_ellipsoids(fits, keep_data = TRUE, add_noise = TRUE, nsim = N)
m2 <- mean_ellipsoids(fits, keep_data = TRUE, add_noise = TRUE, nsim = N^2)

#
plot(els[[1]], col=1, asp=1, xlim=c(-2,2), ylim=c(-2,2))
for(i in 1:5) {
  points(els[[i]], col=i, pch=i)
  plot(fits[[i]], col=i)
}
legend("topleft", legend=1:5, fill=1:5)

plot(m, col=1, lwd=5)
plot(m2, col=2, lwd=5)
#points(m$data, pch=19, cex=0.4)

cat("points :", sapply(els, nrow))
cat("\nestd sd:", round(sqrt(s2<-sapply(lapply(fits, getElement, "ols_fit"), getElement, "s2")),2)  )
cat("\ntrue sd:", round(sds,2))
cat("\n order :", order(s2))
cat("\n m sims:", (m$nsims))
