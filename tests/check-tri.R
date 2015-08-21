tri <- ellipsoid_shape(2, axes = c(1,1,2))
# check coloring
m <- ncol(tri$it)
co <- gray.colors(m)
rescale <- function(v) (v-min(v))/(max(v)-min(v))
co <- gray( rescale(tri$vb[3,c(tri$it)]) )
plot3d(tri, aspect=F, col=co)

w <- triangulation_areas(tri)

co <- rep(gray(rescale(w)), each=3)
plot3d(tri, aspect=F, col=co, lit=F)

