
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
#> reading in input tif: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library/HCR/extdata/2017-01-12_S3A_tsm-Klamath.tif as img_id: 2017-01-12_S3A_tsm-Klamath
#> ...preparing raster to make cost surface
#> ...extracting values from least cost paths
#> ...testing paths for mann kendall
#> connected path to KlamathRiver_25km-ppt11_Klamath River (max: 0.86)
#> .....now saving results in /Users/khondula/Documents/software/HCR/HCRout_klamath_out_2017-01-12_S3A_tsm-Klamath_20260214000219
```

<img src="man/figures/README-example-1.png" width="100%" /><img src="man/figures/README-example-2.png" width="100%" />

    #> img_id: 2017-01-12_S3A_tsm-Klamath
    #> Max HCR point: KlamathRiver_25km-ppt11_Klamath River (0.927)
    #> Closest point: KlamathRiver_25km-ppt4_Damnation Creek

    # most connected site
    klamath_out$which_maxHCR_filtered
    #> [1] "KlamathRiver_25km-ppt11_Klamath River"
