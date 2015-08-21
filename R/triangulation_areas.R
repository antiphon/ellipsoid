#' Compute the areas of polygons per 2D or 3D triangulation
#' 
#' @export
triangulation_areas <- function(tri) {
  if(is(tri, "mesh3d")){
    apply(tri$it, 2, function(ijk){
      vert <- tri$vb[,ijk]
      vert <- t(vert[1:3,])/vert[4,]
      pl <- three_points_to_plane(vert)
      px <- project_to_plane(vert, pl[1:3])
      area_of_2d_polygon(px)
    })
  }
  else{ # 2d triangulation, no need to project
    apply(tri$it, 2, function(ijk){
      vert <- tri$vb[,ijk]
      px <- if(nrow(vert)==3) t(vert[1:2,])/vert[3,] else t(vert)
      area_of_2d_polygon(px)
    })
  }  
}