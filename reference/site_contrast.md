# Environmental contrasts between sites

At-site features: the value of each environmental layer at the source
site minus its value at the destination site. Sites can be points (the
cell value is used) or polygons (the area-weighted mean is used).
Site-level attributes such as meadow area or network degree can be
supplied as a data.frame instead of, or in addition to, rasters.

## Usage

``` r
site_contrast(
  sites,
  pairs,
  env = NULL,
  attributes = NULL,
  id = "site",
  suffix = ".at"
)
```

## Arguments

- sites:

  An `sf` object of sites.

- pairs:

  Data.frame with columns `from` and `to` (directed pairs).

- env:

  Optional `SpatRaster` of environmental layers.

- attributes:

  Optional data.frame with an ID column and one column per site
  attribute.

- id:

  Name of the ID column in `sites` (and `attributes`).

- suffix:

  Appended to the feature names.

## Value

A data.frame with `from`, `to` and one column per feature, each equal to
the source value minus the destination value.

## Examples

``` r
ex <- corridor_example()
pr <- data.frame(from = c("M01", "M12"), to = c("M12", "M01"))
site_contrast(ex$sites, pr, env = ex$env[["snowpack"]])
#>     from  to snowpack.at
#> M01  M01 M12   -496.4841
#> M12  M12 M01    496.4841
```
