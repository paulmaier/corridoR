# Accumulated cost surfaces

The accumulated least cost from each site to every cell of the
landscape. Corridors between two sites are built from the sum of their
surfaces, so this is computed once per site and reused for every pair.

## Usage

``` r
accumulated_cost(tr, sites, id = "site", filename = "")
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

- filename:

  Optional GeoTIFF to write the surfaces to, which keeps memory use low
  for large landscapes.

## Value

A multi-layer `SpatRaster` with one layer per site, named by site ID.
Unreachable cells are `NA`.

## Examples

``` r
# \donttest{
ex <- corridor_example()
tr <- make_transition(ex$resistance)
acc <- accumulated_cost(tr, ex$sites[1:3, ])
# }
```
