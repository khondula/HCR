#' Calculate Hydrologic Connectivity to Rivers
#'
#' @param input_tif filepath for raster map
#' @param sea_point filepath for vector point of interest in coast
#' @param land_points filepath for vector points of river mouths
#' @param land_points_id_col column name for IDs of river mouths
#' @param output_dir folder name to save outputs
#' @param dist_m max distance for river mouths from sea point
#' @param r_min min value for standardizing raster map parameter
#' @param r_max max value for standardizing raster map parameter
#' @param mk_alpha significance level for Mann Kendall test
#' @param mk_tau_threshold tau threshold for Mann Kendall test
#' @param EXP_SCALE_PAR exponential decay parameter for cost surface
#' @param save_outputs save outputs or not
#'
#' @importFrom rlang .data
#' @returns list of results
#' @export
#'
#' @examples
#' klamath_img <- system.file("extdata", "2017-01-12_S3A_tsm-Klamath.tif", package = "HCR")
#' klamath_sea <- system.file("extdata", "sea_pt-Klamath.geojson", package = "HCR")
#' klamath_land <- system.file("extdata", "river_pts-Klamath.geojson", package = "HCR")
#'
#' run_hcr(input_tif = klamath_img,
#'        sea_point = klamath_sea,
#'        land_points = klamath_land,
#'        land_points_id_col = 'pourpoint_id',
#'        output_dir = 'klamath_out',
#'        save_outputs = FALSE)
run_hcr <- function(input_tif,
                    sea_point,
                    land_points,
                    land_points_id_col = 'pourpoint_id',
                    output_dir,
                    dist_m = 25000, # max distance for end points from start pont
                    r_min = 0, # scale for cost surface
                    r_max = 200, # scale for cost surface
                    mk_alpha = 0.05,
                    mk_tau_threshold = 0.5,
                    EXP_SCALE_PAR = 10,
                    save_outputs = TRUE){

  # first make sure all the files exist
  test_file1 <- file.exists(input_tif)
  if(!test_file1){ message('--> input tif does not exist, quitting')
    return(NULL)}

  test_file2 <- file.exists(sea_point)
  if(!test_file2){ message('--> sea_point file does not exist, quitting')
    return(NULL)}

  test_file3 <- file.exists(land_points)
  if(!test_file3){ message('--> land_points file does not exist, quitting')
    return(NULL)}

  # make sure exponential scale parameter is a number
  if(!is.numeric(EXP_SCALE_PAR)){
    EXP_SCALE_PAR <- readr::parse_number(EXP_SCALE_PAR)
  }
  message(glue::glue('EXP_SCALE_PAR is set as {EXP_SCALE_PAR}'))

  #############################
  # READ IN THE PLUME RASTER
  #############################
  img_id <- tools::file_path_sans_ext(basename(input_tif))
  message(glue::glue('reading in input tif: {input_tif} as img_id: {img_id}'))
  input_rast <- terra::rast(input_tif)
  rast_varname <- names(input_rast)
  rast_crs <- terra::crs(input_rast)
  ###########################################################
  # READ IN THE start and end points, convert CRS
  ###########################################################
  start_point_sf <- sf::st_read(sea_point, quiet = TRUE) |>
    sf::st_transform(rast_crs)

  end_points_sf <- sf::st_read(land_points, quiet = TRUE) |>
    sf::st_transform(rast_crs)

  # convert vector data to vect and sp
  start_point_vect <- start_point_sf |> terra::vect()
  start_point_sp <- start_point_sf |> as("Spatial")
  end_points_sp <- end_points_sf |> as("Spatial")

  ###############################################################
  # make sure that the start point is not NA in the plume raster
  ###############################################################

  source_val <- terra::extract(input_rast, start_point_sf)
  test_na <- is.na(source_val[,2])
  if(test_na){ message('--> start point is in NA region of raster, quitting')
    return(NULL)}

  # make output directory
  sub_dir <- glue::glue('HCRout_{output_dir}_{img_id}_{format(Sys.time(), "%Y%M%d%H%M%S")}')
  fs::dir_create(sub_dir)

  ###########################################################
  # FILTER POUR POINTS WITHIN A GIVEN DISTANCE THRESHOLD OF start POINT
  ###########################################################

  # first make a circle around the point
  start_point_25km <- start_point_sf |>
    sf::st_buffer(dist = dist_m) |>
    terra::vect()

  # then filter pour points to that circle
  end_points_25km <- end_points_sf |>
    terra::vect() |>
    terra::intersect(start_point_25km)

  # make an id for these points
  if(is.null(land_points_id_col)){
    end_points_25km$id <- 1:nrow(end_points_25km)
  }
  if(!is.null(land_points_id_col)){
    end_points_25km$id <- end_points_25km[[land_points_id_col]]
  }

  ##############################################
  #### STEP 1) SCALE THE COST SURFACE ##########
  ##############################################

  # EXPONENTIAL SCALING APPROACH TO COST SURFACE

  # Scale input first (0–1)

  # set the min and max values for consistency
  # TSM range: 0 to 200 mg/L

  # # or could get the min and max values from raster
  # r_min <- minmax(tsm_rast)['min',]
  # r_max <- minmax(tsm_rast)['max',]

  r_scaled <- (input_rast - r_min) / (r_max - r_min)

  # Now convert high values to low values
  # using a transformation function (e.g., exponential decay)
  # to Create the resistance surface
  # highest parameter values have lowest resistance
  # default is 10
  # EXP_SCALE_PAR = 20
  # EXP_SCALE_PAR = 50

  my_cost_surface <- exp(-EXP_SCALE_PAR * r_scaled)

  ########################################################
  #### STEP (4) Make transition object from grid values ##
  ########################################################

  # first need to make the raster have only reachable cells

  # drop unreachable area first!!
  # gets rid of any areas on raster that aren't accessible from
  # the source point (e.g. inland bays)
  message('...preparing raster to make cost surface')
  my_cost_surface_reachable <- make_raster_reachable(rast_to_fix = my_cost_surface, my_src_pt = start_point_sp)

  ######################################################
  # then move the watershed pour points onto the raster
  # in the future, could check for maximum distance
  end_pts_moved <- snap_pts_to_rast(pts = end_points_25km, r = my_cost_surface_reachable)
  end_pts_moved <- end_pts_moved |> terra::project(rast_crs)
  ######################################################

  # Find pour point that is the closest for the point of interest
  # calculate eucidlian distances to each point
  eucl_distances <- start_point_vect |>
    terra::distance(end_pts_moved) |>
    as.vector()
  eucl_distances_orig <- start_point_vect |>
    terra::distance(end_points_25km) |>
    as.vector()

  # find which point is the closest (before adjusting to raster)
  min_dist_id <- end_points_25km$id[which.min(eucl_distances_orig)]

  # start making the data frame to save information
  dist_df <- data.frame(pt_id = end_points_25km$id,
                        eucl_dist_orig = eucl_distances_orig) |>
    dplyr::arrange(.data$eucl_dist_orig) |>
    tibble::rowid_to_column(var = 'rank_eucl_dist')

  # make the Transition raster from the cost surface
  # defines the cost of moving between cells
  # easier to move to cells with a lower value
  r1 <- my_cost_surface_reachable |> raster::raster() # convert to old school raster
  tr1 <- gdistance::transition(r1, function(x) 1/mean(x), directions=16)
  tr1 <- gdistance::geoCorrection(tr1, type="c")  # correct for map geometry

  # spatial object type handling
  # make the pour points sf and sp from vect
  end_pts_moved_sf <- end_pts_moved |> sf::st_as_sf()
  end_pts_moved_sp <- end_pts_moved |> sf::st_as_sf() |> as("Spatial")

  ########################################################
  #### STEP (5) Compute least cost paths to each watershed outlet ##
  ########################################################
  # Compute least cost path from start pt (ocean) to each end (pour) point
  # minimum accumulated cost
  # lc_paths is the spatial object of the paths
  # lc_costs is the accumulated costs
  lc_paths <- gdistance::shortestPath(tr1,
                                      origin = start_point_sp,
                                      end_pts_moved_sp,
                                      output="SpatialLines")
  lc_costs <- gdistance::costDistance(tr1,
                                      fromCoords = start_point_sp,
                                      toCoords = end_pts_moved_sp)

  lc_min <- which.min(lc_costs)

  # calculate path lengths
  # convert paths to sf for easy handling
  # copy over ids from the pour pts sp object
  # and get the length of each path
  lc_paths_sf <- sf::st_as_sf(lc_paths) |>
    sf::st_set_crs(terra::crs(input_rast)) |>
    dplyr::mutate(path_id = end_points_25km$id)

  lc_pathlengths <- lc_paths_sf |> sf::st_length()

  df1 <- data.frame(pt_id = end_points_25km$id,
                    path_dist_m = as.numeric(lc_pathlengths),
                    cost = as.vector(lc_costs)) |>
    dplyr::left_join(dist_df, by = c('pt_id')) |>
    dplyr::arrange(.data$cost) |>
    tibble::rowid_to_column(var = 'rank_cost')

  ########################################################
  #### STEP (6) Test paths for positive directional slope ##
  ########################################################
  # then calculate whether there is a monotonic directional gradient
  # towards the shore (for each path)

  message('...extracting values from least cost paths')
  accCost_rast <- gdistance::accCost(tr1, from = start_point_sp)

  accCost_rast <- terra::rast(accCost_rast)
  names(accCost_rast) <- 'acc_cost'
  terra::crs(accCost_rast) <- terra::crs(input_rast)

  lc_paths_list_acc <- 1:nrow(lc_paths_sf) |>
    purrr::map(~sample_lines_and_extract(lc_paths_sf[.x,],
                                         accCost_rast,
                                         d = 150,
                                         path_id_col = 'path_id'))

  lc_paths_list <- 1:nrow(lc_paths_sf) |>
    purrr::map(~sample_lines_and_extract(lc_paths_sf[.x,],
                                         input_rast,
                                         d = 150,
                                         path_id_col = 'path_id'))

  names(lc_paths_list) <- lc_paths_list |>
    purrr::map_chr(~.x$path_id[1])

  vals_end_df <- lc_paths_list |>
    purrr::map(~utils::tail(.x, 1)) |>
    dplyr::bind_rows(.id = 'path_id') |>
    dplyr::select(-.data$x, -.data$y, -.data$dist_along, -.data$seg_dist) |>
    dplyr::rename(end_value = 3,
                  lc_path_dist_m = .data$cum_dist) |>
    dplyr::mutate(pt_id = .data$path_id) |>
    dplyr::select(-.data$path_id)

  message('...testing paths for mann kendall')

  paths_info_list <- lc_paths_list |>
    purrr::map(~test_mannkendall(.x,
                                 yvar = rast_varname,
                                 mk_alpha = mk_alpha,
                                 mk_tau_threshold = mk_tau_threshold))

  paths_summary <- purrr::map(paths_info_list, ~as_tibble(.x[1:2])) |>
    dplyr::bind_rows(.id = 'path_id') |>
    dplyr::mutate(pt_id = .data$path_id) |>
    dplyr::left_join(df1, by = 'pt_id') |>
    dplyr::left_join(vals_end_df, by = 'pt_id') |>
    dplyr::arrange(.data$cost) |>
    dplyr::mutate(cost_div_path = .data$cost/.data$path_dist_m) |>
    dplyr::mutate(HCR = 1 - .data$cost_div_path)

  acc_cost_df <- lc_paths_list_acc |>
    dplyr::bind_rows() |>
    dplyr::select(.data$path_id, .data$x, .data$y, .data$dist_along, .data$cum_dist, .data$acc_cost) |>
    dplyr::mutate(acc_cost = replace(.data$acc_cost, is.infinite(.data$acc_cost), NA))

  lc_paths_df <- lc_paths_list |>
    dplyr::bind_rows() |>
    dplyr::left_join(acc_cost_df, by = c('path_id', 'x', 'y', 'dist_along', 'cum_dist')) |>
    dplyr::mutate(pt_id = .data$path_id) |>
    dplyr::mutate(path_id = as.character(.data$path_id)) |>
    dplyr::left_join(paths_summary, by = c('pt_id', 'path_id'))

  end_points_coords_df <- end_points_25km |>
    sf::st_as_sf() |>
    sf::st_coordinates() |>
    as.data.frame() |>
    dplyr::mutate(pt_id = end_points_25km$id) |>
    dplyr::arrange(.data$Y) |>
    tibble::rowid_to_column(var = 'rank_NorthSouth')

  lc_paths_summary_df <- lc_paths_df |>
    dplyr::select(.data$path_id, .data$pt_id, .data$mk_truepos, .data$mk_tau, .data$rank_cost,
                  .data$end_value, .data$path_dist_m,
                  .data$rank_eucl_dist, .data$HCR) |>
    dplyr::distinct() |>
    dplyr::arrange(-.data$HCR) |>
    tibble::rowid_to_column(var = 'rank_HCR') |>
    dplyr::arrange(.data$HCR) |>
    tibble::rowid_to_column(var = 'rank_revHCR') |>
    dplyr::mutate(mk_connected = .data$mk_truepos & .data$mk_tau >= 0.5) |>
    dplyr::left_join(end_points_coords_df, by = 'pt_id')

  paths_check <- paths_summary |>
    dplyr::filter(.data$mk_truepos) |>
    dplyr::arrange(-.data$HCR)

  if(!nrow(paths_check>0)){
    message(glue::glue('no least cost paths with connected true pos MK > 0.5 (max: {round(max(paths_summary$mk_tau), 2)})'))
  }
  if(nrow(paths_check>0)){
    message(glue::glue('connected path to {paths_check$path_id[1]} (max: {round(max(paths_check$mk_tau[1]), 2)})'))
  }

  maxHCR_hcr <- NA
  which_maxHCR_filtered <- NA

  # if there are points that pass filters!
  n_pathscheck <- nrow(paths_check)

  if(n_pathscheck>0){
    # max HCR
    which_maxHCR_filtered <- paths_check |>
      dplyr::arrange(-.data$HCR) |>
      dplyr::pull(.data$pt_id) |> utils::head(1)

    maxHCR_hcr <- paths_check |>
      dplyr::filter(.data$pt_id == which_maxHCR_filtered) |>
      dplyr::pull(.data$HCR) |> utils::head(1)
  }

  # SAVE things in sub_dir
  # table (all paths data)
  # table (summary by path)
  # points moved
  # paths
  # reachable cost surface

  message(glue::glue('.....now saving results in {sub_dir}'))

  if(save_outputs==TRUE){
    ####################################
    # Make Map Plot
    ####################################
    lc_paths_sf_chk <- lc_paths_sf |>
      dplyr::filter(.data$path_id %in% paths_check$path_id) |>
      dplyr::filter(.data$path_id == which_maxHCR_filtered) |>
      head(1)
    end_pts_connected <- end_pts_moved[which(end_pts_moved$id %in% paths_check$path_id)]
    end_pts_maxHCR <- end_pts_moved[which(end_pts_moved$id %in% which_maxHCR_filtered)]
    end_pts_mindist <- end_pts_moved[which(end_pts_moved$id %in% min_dist_id)]

    mk_tau_df <- lc_paths_df |>
      dplyr::select(.data$path_id, .data$mk_tau) |>
      dplyr::distinct() |>
      dplyr::rename(pourpoint_id = .data$path_id)

    end_pts_moved_vals <- end_pts_moved |>
      terra::values() |>
      dplyr::left_join(mk_tau_df, by = 'pourpoint_id')
    end_pts_moved <- terra::setValues(end_pts_moved, end_pts_moved_vals)

    tau_vals <- end_pts_moved$mk_tau
    my_pal <- colorspace::divergingx_hcl(1000, palette = "Zissou1")
    bg_cols <- scales::col_numeric(
      palette = my_pal,
      domain = c(-1, 1),
      na.color = NA
    )(tau_vals)

    # MAP
    map_filename <- glue::glue('{sub_dir}/plot_map.png')
    png(map_filename, bg = 'white', width = 4, height = 3.5, units = 'in', res = 300)
    raster::plot(my_cost_surface_reachable,
         axes = FALSE,
         box = TRUE,
         colNA = 'gray90',
         mar = c(1,0,3,1),
         cex.main = 0.8,
         main = glue::glue('Max HCR: {which_maxHCR_filtered} ({round(maxHCR_hcr, 2)})'),
         col = scico::scico(100, palette = 'roma', dir = 1))
    terra::lines(lc_paths_sf, col = 'gray')
    terra::points(end_pts_moved, col = 'black', alpha = 0.5, cex = 1, pch = 21, bg = bg_cols)
    terra::lines(lc_paths_sf_chk, col = 'red')
    terra::points(start_point_vect, col = 'red', pch = '*', cex = 2.5)
    terra::points(end_pts_mindist, col = 'black', cex = 0.75, pch = 4, bg = 'black')
    if(nrow(end_pts_maxHCR) > 0){
      terra::points(end_pts_maxHCR, col = 'black', cex = 1, pch = 21, bg = 'red')}
    dev.off()

    ####################################
    # make spaghetti line Plot
    ####################################
    lc_paths_df_maxHCR <- lc_paths_df |>
      dplyr::filter(.data$path_id == which_maxHCR_filtered)
    lc_paths_df_maxHCR_end <- lc_paths_df_maxHCR |>
      dplyr::arrange(-.data$cum_dist) |>
      dplyr::slice_head(n = 1)
    gg1 <- lc_paths_df |>
      ggplot2::ggplot(ggplot2::aes(x = .data$cum_dist/1000, y = .data[[rast_varname]])) +
      ggplot2::geom_line(ggplot2::aes(group = .data$path_id), col = 'gray') +
      ggplot2::geom_line(data = lc_paths_df_maxHCR,
                         ggplot2::aes(group = .data$path_id), col = 'red') +
      ggplot2::geom_point(ggplot2::aes(x = .data$path_dist_m/1000,
                     y = .data$end_value,
                     group = .data$path_id, fill = .data$mk_tau),
                 pch = 21, cex = 3, alpha = 0.9) +
      ggplot2::geom_text(data = lc_paths_df_maxHCR_end,
                         ggplot2::aes(label = glue::glue('HCR: {round(HCR, 2)}   ')),
                size = 3,
                fontface = 'bold',
                hjust = 1) +
      colorspace::scale_fill_continuous_divergingx(name = expression(paste('MK ', tau)), palette = 'Zissou1', mid = 0, limits = c(-1,1)) +
      colorspace::scale_color_continuous_divergingx(palette = 'Zissou1', mid = 0, limits = c(-1,1)) +
      ggplot2::xlab('Least cost path distance (km)') +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    # gg1
    paths_plot_filename <- glue::glue('{sub_dir}/plot_paths.png')
    ggplot2::ggsave(paths_plot_filename, plot = gg1, width = 4, height = 2.5, units = 'in')

  }

  # save things in output directory
  paths_filepath <- glue::glue('{sub_dir}/least_cost_paths.geojson')
  lc_paths_sf |> sf::st_write(paths_filepath)
  end_pts_filepath <- glue::glue('{sub_dir}/land_points_adjusted.geojson')
  end_pts_moved_sf |> sf::st_write(end_pts_filepath)
  lc_paths_df |> readr::write_csv(glue::glue('{sub_dir}/least_cost_paths_data.csv'))
  lc_paths_summary_df |> readr::write_csv(glue::glue('{sub_dir}/least_cost_paths_summary.csv'))


  txt <- glue::glue("img_id: {img_id}
  Max HCR point: {which_maxHCR_filtered} ({round(maxHCR_hcr, 3)})
  Closest point: {min_dist_id}")

  writeLines(txt, glue::glue("{sub_dir}/results_summary.txt"))

  message(txt)

  return_list <- list(img_id = img_id,
                      lc_paths_df = lc_paths_df,
                      lc_paths_summary_df = lc_paths_summary_df,
                      my_cost_surface_reachable = my_cost_surface_reachable,
                      end_pts_moved = end_pts_moved,
                      start_point_vect = start_point_vect,
                      lc_paths_sf = lc_paths_sf,
                      paths_check = paths_check,
                      n_pathscheck = n_pathscheck,
                      which_maxHCR_filtered = which_maxHCR_filtered,
                      min_dist_id = min_dist_id,
                      input_tif = input_tif,
                      sea_point = sea_point,
                      land_points = land_points,
                      land_points_id_col = land_points_id_col)

  }
#' Calculate Hydrologic Connectivity to Rivers
#'
#' @param input_tif filepath for raster map
#' @param sea_point filepath for vector point of interest in coast
#' @param land_points filepath for vector points of river mouths
#' @param land_points_id_col column name for IDs of river mouths
#' @param output_dir folder name to save outputs
#' @param dist_m max distance for river mouths from sea point
#' @param r_min min value for standardizing raster map parameter
#' @param r_max max value for standardizing raster map parameter
#' @param mk_alpha significance level for Mann Kendall test
#' @param mk_tau_threshold tau threshold for Mann Kendall test
#' @param EXP_SCALE_PAR exponential decay parameter for cost surface
#' @param save_outputs save outputs or not
#'
#' @importFrom rlang .data
#' @returns list of results
#' @export
#'
#' @examples
#' klamath_img <- system.file("extdata", "2017-01-12_S3A_tsm-Klamath.tif", package = "HCR")
#' klamath_sea <- system.file("extdata", "sea_pt-Klamath.geojson", package = "HCR")
#' klamath_land <- system.file("extdata", "river_pts-Klamath.geojson", package = "HCR")
#'
#' run_hcr(input_tif = klamath_img,
#'        sea_point = klamath_sea,
#'        land_points = klamath_land,
#'        land_points_id_col = 'pourpoint_id',
#'        output_dir = 'klamath_out',
#'        save_outputs = FALSE)
run_hcr <- function(input_tif,
                    sea_point,
                    land_points,
                    land_points_id_col = 'pourpoint_id',
                    output_dir,
                    dist_m = 25000, # max distance for end points from start pont
                    r_min = 0, # scale for cost surface
                    r_max = 200, # scale for cost surface
                    mk_alpha = 0.05,
                    mk_tau_threshold = 0.5,
                    EXP_SCALE_PAR = 10,
                    save_outputs = TRUE){

  # first make sure all the files exist
  test_file1 <- file.exists(input_tif)
  if(!test_file1){ message('--> input tif does not exist, quitting')
    return(NULL)}

  test_file2 <- file.exists(sea_point)
  if(!test_file2){ message('--> sea_point file does not exist, quitting')
    return(NULL)}

  test_file3 <- file.exists(land_points)
  if(!test_file3){ message('--> land_points file does not exist, quitting')
    return(NULL)}

  # make sure exponential scale parameter is a number
  if(!is.numeric(EXP_SCALE_PAR)){
    EXP_SCALE_PAR <- readr::parse_number(EXP_SCALE_PAR)
  }
  message(glue::glue('EXP_SCALE_PAR is set as {EXP_SCALE_PAR}'))

  #############################
  # READ IN THE PLUME RASTER
  #############################
  img_id <- tools::file_path_sans_ext(basename(input_tif))
  message(glue::glue('reading in input tif: {input_tif} as img_id: {img_id}'))
  input_rast <- terra::rast(input_tif)
  rast_varname <- names(input_rast)
  rast_crs <- terra::crs(input_rast)
  ###########################################################
  # READ IN THE start and end points, convert CRS
  ###########################################################
  start_point_sf <- sf::st_read(sea_point, quiet = TRUE) |>
    sf::st_transform(rast_crs)

  end_points_sf <- sf::st_read(land_points, quiet = TRUE) |>
    sf::st_transform(rast_crs)

  # convert vector data to vect and sp
  start_point_vect <- start_point_sf |> terra::vect()
  start_point_sp <- start_point_sf |> as("Spatial")
  end_points_sp <- end_points_sf |> as("Spatial")

  ###############################################################
  # make sure that the start point is not NA in the plume raster
  ###############################################################

  source_val <- terra::extract(input_rast, start_point_sf)
  test_na <- is.na(source_val[,2])
  if(test_na){ message('--> start point is in NA region of raster, quitting')
    return(NULL)}

  # make output directory
  sub_dir <- glue::glue('{getwd()}/HCRout_{output_dir}_{img_id}_{format(Sys.time(), "%Y%M%d%H%M%S")}')
  fs::dir_create(sub_dir)

  ###########################################################
  # FILTER POUR POINTS WITHIN A GIVEN DISTANCE THRESHOLD OF start POINT
  ###########################################################

  # first make a circle around the point
  start_point_25km <- start_point_sf |>
    sf::st_buffer(dist = dist_m) |>
    terra::vect()

  # then filter pour points to that circle
  end_points_25km <- end_points_sf |>
    terra::vect() |>
    terra::intersect(start_point_25km)

  # make an id for these points
  if(is.null(land_points_id_col)){
    end_points_25km$id <- 1:nrow(end_points_25km)
  }
  if(!is.null(land_points_id_col)){
    end_points_25km$id <- end_points_25km[[land_points_id_col]]
  }

  ##############################################
  #### STEP 1) SCALE THE COST SURFACE ##########
  ##############################################

  # EXPONENTIAL SCALING APPROACH TO COST SURFACE

  # Scale input first (0–1)

  # set the min and max values for consistency
  # TSM range: 0 to 200 mg/L

  # # or could get the min and max values from raster
  # r_min <- minmax(tsm_rast)['min',]
  # r_max <- minmax(tsm_rast)['max',]

  r_scaled <- (input_rast - r_min) / (r_max - r_min)

  # Now convert high values to low values
  # using a transformation function (e.g., exponential decay)
  # to Create the resistance surface
  # highest parameter values have lowest resistance
  # default is 10
  # EXP_SCALE_PAR = 20
  # EXP_SCALE_PAR = 50

  my_cost_surface <- exp(-EXP_SCALE_PAR * r_scaled)

  ########################################################
  #### STEP (4) Make transition object from grid values ##
  ########################################################

  # first need to make the raster have only reachable cells

  # drop unreachable area first!!
  # gets rid of any areas on raster that aren't accessible from
  # the source point (e.g. inland bays)
  message('...preparing raster to make cost surface')
  my_cost_surface_reachable <- make_raster_reachable(rast_to_fix = my_cost_surface, my_src_pt = start_point_sp)

  ######################################################
  # then move the watershed pour points onto the raster
  # in the future, could check for maximum distance
  end_pts_moved <- snap_pts_to_rast(pts = end_points_25km, r = my_cost_surface_reachable)
  end_pts_moved <- end_pts_moved |> terra::project(rast_crs)
  ######################################################

  # Find pour point that is the closest for the point of interest
  # calculate eucidlian distances to each point
  eucl_distances <- start_point_vect |>
    terra::distance(end_pts_moved) |>
    as.vector()
  eucl_distances_orig <- start_point_vect |>
    terra::distance(end_points_25km) |>
    as.vector()

  # find which point is the closest (before adjusting to raster)
  min_dist_id <- end_points_25km$id[which.min(eucl_distances_orig)]

  # start making the data frame to save information
  dist_df <- data.frame(pt_id = end_points_25km$id,
                        eucl_dist_orig = eucl_distances_orig) |>
    dplyr::arrange(.data$eucl_dist_orig) |>
    tibble::rowid_to_column(var = 'rank_eucl_dist')

  # make the Transition raster from the cost surface
  # defines the cost of moving between cells
  # easier to move to cells with a lower value
  r1 <- my_cost_surface_reachable |> raster::raster() # convert to old school raster
  tr1 <- gdistance::transition(r1, function(x) 1/mean(x), directions=16)
  tr1 <- gdistance::geoCorrection(tr1, type="c")  # correct for map geometry

  # spatial object type handling
  # make the pour points sf and sp from vect
  end_pts_moved_sf <- end_pts_moved |> sf::st_as_sf()
  end_pts_moved_sp <- end_pts_moved |> sf::st_as_sf() |> as("Spatial")

  ########################################################
  #### STEP (5) Compute least cost paths to each watershed outlet ##
  ########################################################
  # Compute least cost path from start pt (ocean) to each end (pour) point
  # minimum accumulated cost
  # lc_paths is the spatial object of the paths
  # lc_costs is the accumulated costs
  lc_paths <- gdistance::shortestPath(tr1,
                                      origin = start_point_sp,
                                      end_pts_moved_sp,
                                      output="SpatialLines")
  lc_costs <- gdistance::costDistance(tr1,
                                      fromCoords = start_point_sp,
                                      toCoords = end_pts_moved_sp)

  lc_min <- which.min(lc_costs)

  # calculate path lengths
  # convert paths to sf for easy handling
  # copy over ids from the pour pts sp object
  # and get the length of each path
  lc_paths_sf <- sf::st_as_sf(lc_paths) |>
    sf::st_set_crs(terra::crs(input_rast)) |>
    dplyr::mutate(path_id = end_points_25km$id)

  lc_pathlengths <- lc_paths_sf |> sf::st_length()

  df1 <- data.frame(pt_id = end_points_25km$id,
                    path_dist_m = as.numeric(lc_pathlengths),
                    cost = as.vector(lc_costs)) |>
    dplyr::left_join(dist_df, by = c('pt_id')) |>
    dplyr::arrange(.data$cost) |>
    tibble::rowid_to_column(var = 'rank_cost')

  ########################################################
  #### STEP (6) Test paths for positive directional slope ##
  ########################################################
  # then calculate whether there is a monotonic directional gradient
  # towards the shore (for each path)

  message('...extracting values from least cost paths')
  accCost_rast <- gdistance::accCost(tr1, from = start_point_sp)

  accCost_rast <- terra::rast(accCost_rast)
  names(accCost_rast) <- 'acc_cost'
  terra::crs(accCost_rast) <- terra::crs(input_rast)

  lc_paths_list_acc <- 1:nrow(lc_paths_sf) |>
    purrr::map(~sample_lines_and_extract(lc_paths_sf[.x,],
                                         accCost_rast,
                                         d = 150,
                                         path_id_col = 'path_id'))

  lc_paths_list <- 1:nrow(lc_paths_sf) |>
    purrr::map(~sample_lines_and_extract(lc_paths_sf[.x,],
                                         input_rast,
                                         d = 150,
                                         path_id_col = 'path_id'))

  names(lc_paths_list) <- lc_paths_list |>
    purrr::map_chr(~.x$path_id[1])

  vals_end_df <- lc_paths_list |>
    purrr::map(~utils::tail(.x, 1)) |>
    dplyr::bind_rows(.id = 'path_id') |>
    dplyr::select(-.data$x, -.data$y, -.data$dist_along, -.data$seg_dist) |>
    dplyr::rename(end_value = 3,
                  lc_path_dist_m = .data$cum_dist) |>
    dplyr::mutate(pt_id = .data$path_id) |>
    dplyr::select(-.data$path_id)

  message('...testing paths for mann kendall')

  paths_info_list <- lc_paths_list |>
    purrr::map(~test_mannkendall(.x,
                                 yvar = rast_varname,
                                 mk_alpha = mk_alpha,
                                 mk_tau_threshold = mk_tau_threshold))

  paths_summary <- purrr::map(paths_info_list, ~as_tibble(.x[1:2])) |>
    dplyr::bind_rows(.id = 'path_id') |>
    dplyr::mutate(pt_id = .data$path_id) |>
    dplyr::left_join(df1, by = 'pt_id') |>
    dplyr::left_join(vals_end_df, by = 'pt_id') |>
    dplyr::arrange(.data$cost) |>
    dplyr::mutate(cost_div_path = .data$cost/.data$path_dist_m) |>
    dplyr::mutate(HCR = 1 - .data$cost_div_path)

  acc_cost_df <- lc_paths_list_acc |>
    dplyr::bind_rows() |>
    dplyr::select(.data$path_id, .data$x, .data$y, .data$dist_along, .data$cum_dist, .data$acc_cost) |>
    dplyr::mutate(acc_cost = replace(.data$acc_cost, is.infinite(.data$acc_cost), NA))

  lc_paths_df <- lc_paths_list |>
    dplyr::bind_rows() |>
    dplyr::left_join(acc_cost_df, by = c('path_id', 'x', 'y', 'dist_along', 'cum_dist')) |>
    dplyr::mutate(pt_id = .data$path_id) |>
    dplyr::mutate(path_id = as.character(.data$path_id)) |>
    dplyr::left_join(paths_summary, by = c('pt_id', 'path_id'))

  end_points_coords_df <- end_points_25km |>
    sf::st_as_sf() |>
    sf::st_coordinates() |>
    as.data.frame() |>
    dplyr::mutate(pt_id = end_points_25km$id) |>
    dplyr::arrange(.data$Y) |>
    tibble::rowid_to_column(var = 'rank_NorthSouth')

  lc_paths_summary_df <- lc_paths_df |>
    dplyr::select(.data$path_id, .data$pt_id, .data$mk_truepos, .data$mk_tau, .data$rank_cost,
                  .data$end_value, .data$path_dist_m,
                  .data$rank_eucl_dist, .data$HCR) |>
    dplyr::distinct() |>
    dplyr::arrange(-.data$HCR) |>
    tibble::rowid_to_column(var = 'rank_HCR') |>
    dplyr::arrange(.data$HCR) |>
    tibble::rowid_to_column(var = 'rank_revHCR') |>
    dplyr::mutate(mk_connected = .data$mk_truepos & .data$mk_tau >= 0.5) |>
    dplyr::left_join(end_points_coords_df, by = 'pt_id')

  paths_check <- paths_summary |>
    dplyr::filter(.data$mk_truepos) |>
    dplyr::arrange(-.data$HCR)

  if(!nrow(paths_check>0)){
    message(glue::glue('no least cost paths with connected true pos MK > 0.5 (max: {round(max(paths_summary$mk_tau), 2)})'))
  }
  if(nrow(paths_check>0)){
    message(glue::glue('connected path to {paths_check$path_id[1]} (max: {round(max(paths_check$mk_tau[1]), 2)})'))
  }

  maxHCR_hcr <- NA
  which_maxHCR_filtered <- NA

  # if there are points that pass filters!
  n_pathscheck <- nrow(paths_check)

  if(n_pathscheck>0){
    # max HCR
    which_maxHCR_filtered <- paths_check |>
      dplyr::arrange(-.data$HCR) |>
      dplyr::pull(.data$pt_id) |> utils::head(1)

    maxHCR_hcr <- paths_check |>
      dplyr::filter(.data$pt_id == which_maxHCR_filtered) |>
      dplyr::pull(.data$HCR) |> utils::head(1)
  }

  # SAVE things in sub_dir
  # table (all paths data)
  # table (summary by path)
  # points moved
  # paths
  # reachable cost surface

  message(glue::glue('.....now saving results in {sub_dir}'))

  map_filename <- NA_character_
  paths_plot_filename <- NA_character_
  gg1 <- NULL

  if(save_outputs==TRUE){
    ####################################
    # Make Map Plot
    ####################################
    lc_paths_sf_chk <- lc_paths_sf |>
      dplyr::filter(.data$path_id %in% paths_check$path_id) |>
      dplyr::filter(.data$path_id == which_maxHCR_filtered) |>
      head(1)
    end_pts_connected <- end_pts_moved[which(end_pts_moved$id %in% paths_check$path_id)]
    end_pts_maxHCR <- end_pts_moved[which(end_pts_moved$id %in% which_maxHCR_filtered)]
    end_pts_mindist <- end_pts_moved[which(end_pts_moved$id %in% min_dist_id)]

    mk_tau_df <- lc_paths_df |>
      dplyr::select(.data$path_id, .data$mk_tau) |>
      dplyr::distinct() |>
      dplyr::rename(pourpoint_id = .data$path_id)

    end_pts_moved_vals <- end_pts_moved |>
      terra::values() |>
      dplyr::left_join(mk_tau_df, by = 'pourpoint_id')
    end_pts_moved <- terra::setValues(end_pts_moved, end_pts_moved_vals)

    tau_vals <- end_pts_moved$mk_tau
    my_pal <- colorspace::divergingx_hcl(1000, palette = "Zissou1")
    bg_cols <- scales::col_numeric(
      palette = my_pal,
      domain = c(-1, 1),
      na.color = NA
    )(tau_vals)

    # MAP
    draw_hcr_map <- function() {
      raster::plot(my_cost_surface_reachable,
                   axes = FALSE,
                   box = TRUE,
                   colNA = 'gray90',
                   mar = c(1,0,3,1),
                   cex.main = 0.8,
                   main = glue::glue('Max HCR: {which_maxHCR_filtered} ({round(maxHCR_hcr, 2)})'),
                   col = scico::scico(100, palette = 'roma', dir = 1))
      terra::lines(lc_paths_sf, col = 'gray')
      terra::points(end_pts_moved, col = 'black', alpha = 0.5, cex = 1, pch = 21, bg = bg_cols)
      terra::lines(lc_paths_sf_chk, col = 'red')
      terra::points(start_point_vect, col = 'red', pch = '*', cex = 2.5)
      terra::points(end_pts_mindist, col = 'black', cex = 0.75, pch = 4, bg = 'black')
      if(nrow(end_pts_maxHCR) > 0){
        terra::points(end_pts_maxHCR, col = 'black', cex = 1, pch = 21, bg = 'red')}
    }

    draw_hcr_map()

    map_filename <- glue::glue('{sub_dir}/plot_map.png')

    png(map_filename, bg = 'white', width = 4, height = 3.5, units = 'in', res = 300)
    draw_hcr_map()
    dev.off()

    ####################################
    # make spaghetti line Plot
    ####################################
    lc_paths_df_maxHCR <- lc_paths_df |>
      dplyr::filter(.data$path_id == which_maxHCR_filtered)
    lc_paths_df_maxHCR_end <- lc_paths_df_maxHCR |>
      dplyr::arrange(-.data$cum_dist) |>
      dplyr::slice_head(n = 1)
    gg1 <- lc_paths_df |>
      ggplot2::ggplot(ggplot2::aes(x = .data$cum_dist/1000, y = .data[[rast_varname]])) +
      ggplot2::geom_line(ggplot2::aes(group = .data$path_id), col = 'gray') +
      ggplot2::geom_line(data = lc_paths_df_maxHCR,
                         ggplot2::aes(group = .data$path_id), col = 'red') +
      ggplot2::geom_point(ggplot2::aes(x = .data$path_dist_m/1000,
                                       y = .data$end_value,
                                       group = .data$path_id, fill = .data$mk_tau),
                          pch = 21, cex = 3, alpha = 0.9) +
      ggplot2::geom_text(data = lc_paths_df_maxHCR_end,
                         ggplot2::aes(label = glue::glue('HCR: {round(HCR, 2)}   ')),
                         size = 3,
                         fontface = 'bold',
                         hjust = 1) +
      colorspace::scale_fill_continuous_divergingx(name = expression(paste('MK ', tau)), palette = 'Zissou1', mid = 0, limits = c(-1,1)) +
      colorspace::scale_color_continuous_divergingx(palette = 'Zissou1', mid = 0, limits = c(-1,1)) +
      ggplot2::xlab('Least cost path distance (km)') +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    # gg1
    paths_plot_filename <- glue::glue('{sub_dir}/plot_paths.png')
    ggplot2::ggsave(paths_plot_filename, plot = gg1, width = 4, height = 2.5, units = 'in')

    print(gg1)

  }

  # save things in output directory
  paths_filepath <- glue::glue('{sub_dir}/least_cost_paths.geojson')
  lc_paths_sf |> sf::st_write(paths_filepath)
  end_pts_filepath <- glue::glue('{sub_dir}/land_points_adjusted.geojson')
  end_pts_moved_sf |> sf::st_write(end_pts_filepath)
  lc_paths_df |> readr::write_csv(glue::glue('{sub_dir}/least_cost_paths_data.csv'))
  lc_paths_summary_df |> readr::write_csv(glue::glue('{sub_dir}/least_cost_paths_summary.csv'))


  txt <- glue::glue("img_id: {img_id}
  Max HCR point: {which_maxHCR_filtered} ({round(maxHCR_hcr, 3)})
  Closest point: {min_dist_id}")

  writeLines(txt, glue::glue("{sub_dir}/results_summary.txt"))

  message(txt)

  return_list <- list(img_id = img_id,
                      map_file = map_filename,
                      paths_plot = gg1,
                      lc_paths_df = lc_paths_df,
                      lc_paths_summary_df = lc_paths_summary_df,
                      my_cost_surface_reachable = my_cost_surface_reachable,
                      end_pts_moved = end_pts_moved,
                      start_point_vect = start_point_vect,
                      lc_paths_sf = lc_paths_sf,
                      paths_check = paths_check,
                      n_pathscheck = n_pathscheck,
                      which_maxHCR_filtered = which_maxHCR_filtered,
                      min_dist_id = min_dist_id,
                      input_tif = input_tif,
                      sea_point = sea_point,
                      land_points = land_points,
                      land_points_id_col = land_points_id_col)

}




