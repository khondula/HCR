#' Make raster reachable
#'
#' @param rast_to_fix A raster with isolated water pixels
#' @param my_src_pt Point vector object for ocean point of interest
#' @param src_pt_index Row of ocean point of interest in my_src_pt, default = 1
#'
#' @returns Raster object with pixels removed that arent connected to the source point
#' @export
#'
#' @examples
#' # reachable_raster <- make_raster_reachable(rat)
make_raster_reachable <- function(rast_to_fix, my_src_pt, src_pt_index = 1){
  src_pt <- my_src_pt[src_pt_index,]
  r1 <- rast_to_fix |> raster::raster() # convert to old school raster
  tr1 <- gdistance::transition(r1, function(x) 1/mean(x), directions=16)
  tr1 <- gdistance::geoCorrection(tr1, type="c")  # correct for map geometry

  # 1) accumulated least-cost distance from source
  acc <- gdistance::accCost(tr1, src_pt)   # acc is a RasterLayer of accumulated cost (Inf = unreachable)

  # 2) create binary reachable mask: 1 = reachable, NA = unreachable
  reachable <- acc
  reachable[is.infinite(reachable[])] <- NA    # make unreachable NA
  reachable[!is.na(reachable[])] <- 1         # make reachable 1s

  # 3) label connected components (8-neighbourhood)
  cl <- raster::clump(reachable, directions = 8)

  # 4) find which component contains the source
  src_cell <- raster::cellFromXY(cl, raster::coordinates(src_pt))
  src_comp <- cl[src_cell]
  if (is.na(src_comp)) stop("Source cell is not in any reachable component (it may be isolated).")

  # 5) build a mask that keeps only that component
  keep_mask <- cl
  keep_mask[keep_mask != src_comp] <- NA
  keep_mask[keep_mask == src_comp] <- 1

  # 6) apply mask to your original cost raster (unreachable -> NA)
  cost_reachable <- raster::mask(r1, keep_mask)
  cost_reachable_rast <- terra::rast(cost_reachable)
  return(cost_reachable_rast)

}

#' Snap points to raster
#'
#' @param pts vector points of watershed pour points
#' @param r raster of cost surface
#'
#' @returns vector points moved onto valid pixels of raster
#' @export
#'
#' @examples
#' #   end_pts_moved <- snap_pts_to_rast(pts = end_points_25km, r = my_cost_surface_reachable)
snap_pts_to_rast <- function(pts, r){
  # Get coordinates of non-NA cells
  valid_cells <- which(!is.na(raster::values(r)))
  xy_valid <- raster::xyFromCell(r, valid_cells)

  # Extract coordinates of your points
  xy_pts <- raster::geom(pts)[, c("x", "y")]
  # xy_pts <- geom(pts)[24, c("x", "y"), drop = FALSE] # to test for one point
  #  Find nearest valid raster cell for each point
  nn <- FNN::get.knnx(data = xy_valid, query = xy_pts, k = 1)

  #  Get new coordinates (nearest valid cell centers)
  xy_new <- xy_valid[nn$nn.index[, 1], ]
  colnames(xy_new) <- c("x", "y")

  # Combine new coordinates with original attributes
  df_new <- cbind(xy_new, as.data.frame(pts))

  # Create new point vector with same CRS and attributes
  pts_moved <- terra::vect(df_new, geom = c("x", "y"), crs = terra::crs(r))

  return(pts_moved)
}

