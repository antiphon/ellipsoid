#' test OLS

library(devtools)
load_all(".")

axes <- c(1,2)
R <- sphere::rotationMatrix(az=0.1)[-3,-3]
x <- rellipsoid(50, axes, noise.sd = 0.1, R = R, center = c(2,0))

e <- ellipsoid_OLS(x, origin=F)

plot(x, asp=1)
#plot.ellipsoid(e)

plot(e, add = T)
