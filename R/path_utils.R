#' Sample lines and extract
#'
#' @param lines_sf line vector to extract along
#' @param rast raster to extract from
#' @param d distance to sample line along
#' @param path_id_col column name for path ID
#' @param drop_na should NA values be dropped
#'
#' @returns table with extractex values
#' @export
#'
#' @examples
#' # sample_lines_and_extract(lc_paths_sf, input_rast, d = 150, path_id_col = 'path_id')
sample_lines_and_extract <- function(lines_sf, rast, d = 150, path_id_col = NULL, drop_na = TRUE) {
  # helper function: sample a single line at regular metre spacing d (metres)
  # ensure raster CRS string (terra) and transform lines to that CRS
  rast_crs <- terra::crs(rast)            # terra CRS string
  lines_proj <- sf::st_transform(lines_sf, crs = rast_crs)

  # helper: sample a single line at metre spacing d and extract raster values
  sample_one <- function(line_sf, pid) {
    L <- as.numeric(sf::st_length(line_sf))
    if (L == 0) return(tibble::tibble())   # degenerate

    # metre positions including endpoint
    pos_m <- seq(0, L, by = d)
    if (utils::tail(pos_m, 1) < L) pos_m <- c(pos_m, L)
    frac <- pos_m / L

    # sample points (projected CRS)
    pts_sfc <- sf::st_line_sample(line_sf, sample = frac)
    # cast MULTIPOINT -> individual POINT geometries (one per sampled position)
    pts_pts <- sf::st_cast(pts_sfc, "POINT")
    coords <- sf::st_coordinates(pts_pts)

    # build tibble of coordinates + distance along
    df_pts <- tibble::tibble(
      path_id = pid,
      x = coords[, "X"],
      y = coords[, "Y"],
      dist_along = pos_m
    )
    # Euclidean seg / cumulative distances (metres)
    dx <- c(0, diff(df_pts$x))
    dy <- c(0, diff(df_pts$y))
    df_pts <- df_pts |>
      dplyr::mutate(seg_dist = sqrt(dx^2 + dy^2),
             cum_dist = cumsum(.data$seg_dist))

    # convert to terra SpatVector for extraction
    pts_vect <- terra::vect(df_pts,
                     geom = c("x", "y"),
                     crs = rast_crs)
    vals <- terra::extract(rast, pts_vect)
    # terra::extract returns a data.frame; first column often "ID" (point index)
    # drop ID to avoid duplication, then bind with df_pts
    if ("ID" %in% names(vals)) vals$ID <- NULL
    # ensure row-order matches (terra returns same order as input)
    out <- dplyr::bind_cols(df_pts, tibble::as_tibble(vals))

    if (drop_na) {
      raster_cols <- setdiff(names(out), names(df_pts))
      out <- out |> dplyr::filter(dplyr::if_all(tidyselect::all_of(raster_cols), ~ !is.na(.x)))
    }

    out
  }

  # determine path ids (preserve user-specified column if provided)
  if (!is.null(path_id_col) && path_id_col %in% names(lines_proj)) {
    ids <- lines_proj[[path_id_col]]
  } else {
    ids <- seq_len(nrow(lines_proj))
  }

  # iterate over features
  res_list <- purrr::map2(lines_proj |> split(seq_len(nrow(lines_proj))),
                   ids,
                   ~ sample_one(.x, .y))

  # combine into single tibble (preserve order)
  res_df <- dplyr::bind_rows(res_list)
  return(res_df)
}
