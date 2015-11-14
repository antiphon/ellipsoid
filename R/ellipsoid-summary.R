#############################################################
#' Summarise an ellipsoid
#' 
#' @export
summary.ellipsoid <- function(x, ...){
  print(x)
  ang <- round(x$rot_angle, 3)
  angd <- round(x$rot_angle * 180/pi, 1)
  semi_axes <- x$semi_axes
  semi_axes_rel <- round(semi_axes/semi_axes[1], 3)
  angtxt <- ifelse(x$dim==2, format(ang), paste(c("heading", "attitude", "bank"), format(ang), collapse=" "))
  angtxtd <- ifelse(x$dim==2, format(angd), paste(c("heading", "attitude", "bank"), format(angd), collapse=" "))
  cat("\nEstimates:\n Center:\t \t ", paste0("(", paste0(format(x$center), collapse=", "), ")\n"))
  cat(" Semi-axes lengths (absolute):\t ", paste0(format(semi_axes), collapse=" : "), "\n")
  cat(" Semi-axes lengths (relative):\t ", paste0(format(semi_axes_rel), collapse=" : "), "\n")
  
  cat(" Rotation angles (rad):\t ", angtxt,"\n")
  cat(" Rotation angles (deg):\t ", angtxtd,"\n")
  cat(" Error variance: \t ", x$ols_fit$s2, "\n")
  invisible(NULL)
}
