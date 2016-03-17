# test the 2d plot of 3D ellipse

library(devtools)
load_all(".")

rr <- rellipsoid(1000, c(1,2,1), 
                 R=sphere::rotationMatrix(ax=runif(1,-pi/3,pi/5)))

e <- ellipsoid_OLS(rr)

par(mfrow=c(3,1))

for(i in 1:3) plot(e, i=i, add=F, levels=12)

library(rgl)
plot3d(cbind(x=0,y=0,z=0))
plot(e, alpha=.5)
