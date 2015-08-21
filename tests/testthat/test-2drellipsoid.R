# test for 2D rellipsoid

test_that("2D ellipsoids sampling test, roughly uniform",
{
  axes <- c(1,2)
  set.seed(1)
  x <- rellipsoid(10000, axes=axes)
  pieces <- 100
  theta <- seq(0, 2 * pi, length=pieces)
  xy <- cbind(axes[1] * cos(theta), axes[2] * sin(theta))
  # the surface integrals ~ distance between consecutive points
  integ <- rowSums((xy[-1,] - xy[-pieces,])^2)
  den <- integ/sum(integ)
  # bin
  thet <- atan2(x[,2],x[,1]) # angles of data
  thet[thet<0] <- 2*pi + thet[thet<0]
  z <- cut(thet, theta, labels=FALSE)
  ncount <- tabulate(z, nbins = pieces-1)
  f <- 1/(ncount/sum(ncount))
  actual <- round(sqrt(sum(f-den)))
  # the number of points should inversely related to the size of the surface
  expected <- 111
  expect_equal(actual, expected)
})
