# Summarize environmental layers within corridors

For every pair of sites, summarizes each environmental layer inside the
pair's corridor. With least cost corridors (`acc` and `q`), each cell
counts in proportion to its corridor weight, so likely routes matter
more than marginal ones. With path buffers, cells are weighted by the
share of the cell inside the buffer.

## Usage

``` r
corridor_extract(
  env,
  pairs,
  acc = NULL,
  q = 0.05,
  buffers = NULL,
  stat = "mean",
  progress = interactive()
)
```

## Arguments

- env:

  A `SpatRaster` of environmental layers. Layers are resampled to the
  corridor grid as needed (bilinear for continuous layers).

- pairs:

  Data.frame with columns `from` and `to`.

- acc, q:

  Accumulated cost surfaces from
  [`accumulated_cost()`](https://paulmaier.github.io/corridoR/reference/accumulated_cost.md)
  and the corridor threshold for
  [`lcc_weights()`](https://paulmaier.github.io/corridoR/reference/lcc_weights.md).
  Use these for least cost corridors.

- buffers:

  Alternatively, an `sf` layer of buffered paths from
  [`path_buffers()`](https://paulmaier.github.io/corridoR/reference/path_buffers.md)
  with columns `from` and `to`.

- stat:

  How to summarize each layer: `"mean"` (weighted mean, the default) or
  `"sum"` (weighted sum, for counts such as stream or trail crossings).
  A named vector sets the statistic per layer.

- progress:

  Show a progress bar.

## Value

A data.frame with `from`, `to` and one column per layer of `env`.

## Examples

``` r
# \donttest{
ex <- corridor_example()
tr <- make_transition(ex$resistance)
pr <- data.frame(from = c("M01", "M02"), to = c("M12", "M05"))
acc <- accumulated_cost(tr, ex$sites[ex$sites$site %in% unlist(pr), ])
corridor_extract(ex$env, pr, acc = acc, q = 0.05)
#>   from  to snowpack   runoff moisture summer_temp    slope    forest
#> 1  M01 M12 702.5725 474.9524 7.035619    22.80098 3.897698 0.9004108
#> 2  M02 M05 720.7984 483.6719 7.086381    22.82677 3.646871 0.9465355
# }
```
