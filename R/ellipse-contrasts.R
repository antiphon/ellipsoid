#' Testing contrast for ellipses
#' 
#' Give this to confint.ellipsoid to compute the contrast for
#'  testing semi_a - semi_b = 0
#' @param out do outlier removal, the quadratic simulations is unstable.
#' @param range range of outlier removal. See boxplot for details of the formula.
#' 
#' 
#' @export

ellipse_contrast_2d <- function(els, out=TRUE, range=1.85){
  
  if(!ncol(els[[1]]$A)==2) stop("Only for 2D.")
  
  R <- lapply(lapply(els, getElement, "A"), ellipse_solve_rota)
  axes <- t(sapply(R, function(b) b$axes ))
  angs <- sapply(lapply(R, getElement, "R"), function(R) {f <- R %*% c(1,0)
  atan2(f[2],f[1])})
  angs[angs <0 ] <- angs[ angs<0 ] + pi
  ab <- sapply(1:length(angs), function(i) 
    if(angs[i] < 3*pi/4 & angs[i] > pi/4) axes[i,1:2] else axes[i,2:1]  )
  
  e<- ab[1,]-ab[2,]
  e <- e[!is.na(e) & !is.infinite(e)]
  # check outliers
  if(out){
    m <- median(e, na.rm=TRUE)
    out <- abs(e-m) > range * abs(quantile(e-m, prob=3/4))
    e<-e[!out]
  }
  e
  
}

