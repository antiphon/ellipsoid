#' test ALS

library(devtools)
load_all(".")

axes <- c(1, 2, 3)
R <- sphere::rotationMatrix(az=pi/4)
x <- rellipsoid(250, axes, noise.sd = 0.1, R = R)

ee <- ellipsoid_OLS(x)
se <- ee$ols_fit$s2
e <- ellipsoid_ALS(x, s2=seq(se*.1,se, length=20))

#plot(x, asp=1)
#plot.ellipsoid(e)

summary(e)
cat(c(pi/4, pi-pi/4))
#plot(e)

par(mfrow=c(2,3))
persp(ee, triptych = T)
persp(e, triptych = T)
