library(testthat)
#library(devtools)
library(ellipsoid)

options(testthat.use_colours = FALSE)
test_check("ellipsoid")

# 
# with_envvar(
#   c(LANG = "en_US"), test_package("ellipsoid")
#)