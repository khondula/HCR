
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HCR

<!-- badges: start -->

<!-- badges: end -->

Code for implementing hydrologic connectivity to rivers.

## Installation

You can install the development version of HCR from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("khondula/HCR")
```

## Example

This is a basic example to calculate HCR:

``` r
library(HCR)
klamath_img_path <- system.file("extdata", "2017-01-12_S3A_tsm-Klamath.tif", package = "HCR")
klamath_sea_path <- system.file("extdata", "sea_pt-Klamath.geojson", package = "HCR")
klamath_land_path <- system.file("extdata", "river_pts-Klamath.geojson", package = "HCR")


klamath_out <- run_hcr(input_tif = klamath_img_path,
        sea_point = klamath_sea_path,
        land_points = klamath_land_path,
        land_points_id_col = 'pourpoint_id',
        output_dir = 'klamath_out')
#> EXP_SCALE_PAR is set as 10
#> reading in input tif: C:/Users/khondula/AppData/Local/R/win-library/4.5/HCR/extdata/2017-01-12_S3A_tsm-Klamath.tif as img_id: 2017-01-12_S3A_tsm-Klamath
#> ...preparing raster to make cost surface
#> ...extracting values from least cost paths
#> ...testing paths for mann kendall
#> connected path to KlamathRiver_25km-ppt11_Klamath River (max: 0.86)
#> .....now saving results in HCRout_klamath_out_2017-01-12_S3A_tsm-Klamath_20265514025515
#> Writing layer `least_cost_paths' to data source 
#>   `HCRout_klamath_out_2017-01-12_S3A_tsm-Klamath_20265514025515/least_cost_paths.geojson' using driver `GeoJSON'
#> Writing 8 features with 1 fields and geometry type Line String.
#> Writing layer `land_points_adjusted' to data source 
#>   `HCRout_klamath_out_2017-01-12_S3A_tsm-Klamath_20265514025515/land_points_adjusted.geojson' using driver `GeoJSON'
#> Writing 8 features with 17 fields and geometry type Point.

# most connected site
klamath_out$which_maxHCR_filtered
#> [1] "KlamathRiver_25km-ppt11_Klamath River"
```
