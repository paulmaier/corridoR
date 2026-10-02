# Least cost paths between sites

Least cost paths between sites

## Usage

``` r
least_cost_paths(
  tr,
  sites,
  id = "site",
  pairs = NULL,
  dem = NULL,
  max_length = Inf
)
```

## Arguments

- tr:

  A transition layer from
  [`make_transition()`](https://paulmaier.github.io/corridoR/reference/make_transition.md).

- sites:

  An `sf` object of points or polygons (polygon centroids are used) with
  an ID column.

- id:

  Name of the ID column in `sites`.

- pairs:

  Optional data.frame with columns `from` and `to`. By default every
  unordered pair of sites.

- dem:

  Optional elevation `SpatRaster`. If given, path length is measured
  over the terrain surface (3D), as in Maier et al. (2022).

- max_length:

  Drop paths longer than this (meters).

## Value

An `sf` object of lines with columns `from`, `to`, `euclidean` and
`length`.

## Examples

``` r
ex <- corridor_example()
tr <- make_transition(ex$resistance)
paths <- least_cost_paths(tr, ex$sites, pairs = data.frame(from = "M01", to = "M12"),
                          dem = ex$dem)
paths
#> Simple feature collection with 1 feature and 4 fields
#> Geometry type: LINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: 3500 ymin: 3100 xmax: 30100 ymax: 23700
#> Projected CRS: +proj=tmerc +lat_0=0 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#>   from  to euclidean                       geometry   length
#> 1  M01 M12   19636.7 LINESTRING (3500 19500, 370... 57765.03
```
