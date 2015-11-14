#' test confint

library(devtools)
load_all(".")

axes <- c(1,2)
R <- sphere::rotationMatrix(az=0.1)[-3,-3]

xx <- rellipsoid(100, axes, noise.sd = 0.2, R = R)

x <- ellipsoid_OLS(xx)



t0 <- system.time( ci <- confint.ellipsoid_old(x) )
t1 <- system.time( ci2 <- confint(x) )

print(rbind(t0,t1))

par(mfrow=c(1,2))
plot(ci)
plot(ci2)
