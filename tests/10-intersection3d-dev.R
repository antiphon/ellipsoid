# test ellipsoid intersection with plane
library(rgl)
library(devtools)
load_all(".")


R <- sphere::rotationMatrix(.2,.1,0)
sa <- c(1, 2.5, 5)/2
x <- as_ellipsoid(sa, R, c(0,0,0))
# normal
#n <- c(0,1,0); n <- n/sqrt(sum(n^2))
#q <- c(0,0,0)
#p <- intersect_ellipsoid_plane(x, n, q)
# check
#basis <- p$basis
#print(t(basis)%*%basis)
###########################################
# for one of the basis
di <- 1
n <- c(0, 0, 0); n[di] <- 1
q <- c(0, 0, 0) + x$center
p <- intersect_ellipsoid_plane(x, n, q)
#########################################
#a <- seq(0, 2*pi, l=100)
#xy <- cbind(cos(a), sin(a)) %*% t(R%*%diag(p$semi_axes))
#########################################
if(1){
plot3d(cbind(x=0, y=0, z=0))
plot(x, alpha=.5, col=7, N=3)
A <- p$semi_axes[1]
B <- p$semi_axes[2]
r <- p$basis[,1]
s <- p$basis[,2]
m <- p$center3d
gx <-seq(-2, 2, l=150)
g <- expand.grid(gx,gx)
xx <- t( apply(g,1, function(gi) m + gi[1]*r+ gi[2]*s)  )
points3d(xx, alpha=.5, size=1)
lines3d(rbind(m,m+n), col=2)
lines3d(rbind(m,m+r*A), col=3)
lines3d(rbind(m,m+s*B), col=4)
points3d(m)
# 
plot(x, i=di, col=3, add=F, levels=10)

}






# rubbish

# summary(el)
# R <- el$rot
# D <- diag(1/el$semi_axes)
# # orthongonal vectors.
# n0 <- c(0,0,1); n0 <- n0/sqrt(sum(n0^2))
# q <- c(c(0,0,.3)%*%R)
# n <- c(n0%*%R) 
# # random unit
# u <- sphere::runifsphere(1)#c(1,0,0)%*%R
# r <- cross(u, n)
# s <- cross(n, r)
# # helper
# Df <- function(a,b) t(D%*%a)%*%(D%*%b)
# dot <- function(a,b) t(a)%*%b
# # reorient
# Drs <-Df(r,s)
# Drr <-Df(r,r)
# Dss <-Df(s,s)
# om <- 0.5 * atan(2*(Drs)/(Drr-Dss))
# r <- cos(om)*r + sin(om)*s
# s <- -sin(om)*r + cos(om)*s
# # checks, should be 0
# print(c(Df(r,s)))
# print(c(dot(n,r), dot(n,s), dot(s,r)))
# # 
# Dqr <- Df(q,r)
# Dqs <- Df(q,s)
# Drr <- Df(r,r)
# Dss <- Df(s,s)
# cent <- c( -Dqr/Drr, -Dqs/Dss)
# d <- Df(q,q) - Dqr^2/Drr-Dqs^2/Dss
# semi_axes <- sqrt(c(a= (1-d)/Drr, b=(1-d)/Dss ))
# illustrate
# a <- seq(0,pi*2,l=100)
# xy0 <- cbind(cos(a), sin(a)) %*% diag(semi_axes)
# xy <- xy0
# plot(xy, asp=1, type="l")
# plot(x, N=3, alpha=.2, col=4)
# # plot in 3d
# m <- q + cent[1]*r + cent[2]*s
# xyz0 <- cbind(xy, 0)
# A <- semi_axes[1]
# B <- semi_axes[2]
# xyz <- t(sapply(a, function(ai) m + A*cos(ai)*r+B*sin(ai)*s  ))
# Ri <- t(R)
# xyz <- xyz %*% Ri
# 
# 
# #plot3d(er, aspect=F, alpha=0.3)
# 
# plot(el, N=3, alpha=.2)
# 
# lines3d(xyz, col=4)
# print(det(R))
# lines3d(rbind(q,q+n)%*%Ri, col=1)
# lines3d(rbind(q,q+s)%*%Ri, col=2)
# lines3d(rbind(q,q+r)%*%Ri, col=3)
# 
