#' test ALS 2d

library(devtools)
load_all(".")

axes <- c(1, 5)
R <- sphere::rotationMatrix(az=pi/4)[-3,-3]
x <- rellipsoid(50, axes, noise.sd = s<-0.2, R = R)

ee <- ellipsoid_OLS(x)
se <- ee$ols_fit$s2
e <- ellipsoid_ALS(x, s2=seq(se*.1, se, length=10))

summary(ee)
summary(e)
cat("Noise s2:", s^2,"\nRot was:", c(pi/4, pi-pi/4),"\n")

print(ee$ols_fit$varcov)
print(e$ols_fit$varcov)

confint(ee)
confint(e)

par(mfrow=c(1,1))
plot(x, asp=1, xlim=c(-5,5)*2)
plot(e)
plot(ee, col=2)
