#' test OLS

library(devtools)
load_all(".")

axes <- c(1, 2, 3)
R <- sphere::rotationMatrix(az=0.5)
x <- rellipsoid(50, axes, noise.sd = 0.1, R = R)

e <- ellipsoid_OLS(x, origin=T)

#plot(x, asp=1)
#plot.ellipsoid(e)

summary(e)

#plot(e)

par(mfrow=c(3,1))
persp(e, triptych = T)
