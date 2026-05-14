
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HCR

<!-- badges: start -->

<!-- badges: end -->

Code for implementing hydrologic connectivity to rivers.

## Installation

You can install the development version of HCR like so:

``` r
library(devtools)
devtools::install_github('khondula/HCR')
```

## Example

This is a basic example to calculate HCR:

``` r
library(HCR)
klamath_img_path <- system.file("extdata", "2017-01-12_S3A_tsm-Klamath.tif", package = "HCR")
klamath_sea_path <- system.file("extdata", "sea_pt-Klamath.geojson", package = "HCR")
klamath_land_path <- system.file("extdata", "river_pts-Klamath.geojson", package = "HCR")


run_hcr(input_tif = klamath_img_path,
        sea_point = klamath_sea_path,
        land_points = klamath_land_path,
        land_points_id_col = 'pourpoint_id',
        output_dir = 'klamath_out')
#> --> input tif does not exist, quitting
#> NULL
```
