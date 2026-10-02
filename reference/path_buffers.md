# Buffers around least cost paths

Simple corridor bandwidths: a fixed distance on either side of each
least cost path.

## Usage

``` r
path_buffers(paths, width = c(100, 500))
```

## Arguments

- paths:

  Least cost paths from
  [`least_cost_paths()`](https://paulmaier.github.io/corridoR/reference/least_cost_paths.md).

- width:

  Buffer width(s) in meters.

## Value

A named list of `sf` polygon layers, one per width.

## Examples

``` r
# \donttest{
ex <- corridor_example()
tr <- make_transition(ex$resistance)
p <- least_cost_paths(tr, ex$sites, pairs = data.frame(from = "M01", to = "M12"))
b <- path_buffers(p, c(100, 500))
# }
```
