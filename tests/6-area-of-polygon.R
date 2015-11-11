# test area of 2d polygon

library(devtools)
load_all(".")

library(spatstat)
x <- runifpoint(100)
e <- ripras(x)
p <- with(e$bdry[[1]], cbind(x,y))

polygon(p)

a <- c(area(e), area_of_2d_polygon(p))
print(a)
