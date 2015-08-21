#' Compute the area of 2D polygon given by ordered border points
#' 
#' So called shoelace method.
#' 
#' @export
area_of_2d_polygon <- function(px){
  px <- rbind(px,px[nrow(px),])
  n <- nrow(px)
  0.5 * abs( sum(px[-n,1]*px[-1,2]) - sum(px[-1,1]*px[-n,2]) )
}


