# Arrows for a directional map

Summarizes a map from
[`map_pairwise()`](https://paulmaier.github.io/corridoR/reference/map_pairwise.md)
on a coarser grid, for drawing arrows. Direction follows the summed
vectors. By default length does too, so arrows are long where corridors
agree on a direction and shrink where opposing flows cancel;
`length_by = "magnitude"` scales them by the mean size of change
instead.

## Usage

``` r
shift_arrows(
  map,
  cells = 25,
  length = NULL,
  length_by = c("net", "magnitude"),
  scale = c("linear", "sqrt"),
  min_quantile = 0,
  min_length = 0
)
```

## Arguments

- map:

  A directional map from
  [`map_pairwise()`](https://paulmaier.github.io/corridoR/reference/map_pairwise.md).

- cells:

  Grid size, in cells of `map`.

- length:

  Length of the longest arrow, in map units. Defaults to about 80% of a
  grid cell.

- length_by:

  `"net"` (length of the summed vector) or `"magnitude"` (mean size of
  change in the grid cell).

- scale:

  `"linear"` lengths, or `"sqrt"` to compress the range so a few very
  large changes do not shrink every other arrow to a dot.

- min_quantile:

  Drop arrows in the weakest cells (below this quantile of magnitude).

- min_length:

  Shortest arrow, as a share of `length`. With the default of 0, arrow
  length is proportional to magnitude, so arrows shrink to nothing where
  nothing changes.

## Value

A data.frame with `x`, `y`, `xend`, `yend` and `magnitude`, ready for
[`ggplot2::geom_segment()`](https://ggplot2.tidyverse.org/reference/geom_segment.html)
or [`graphics::arrows()`](https://rdrr.io/r/graphics/arrows.html).

## Examples

``` r
# \donttest{
ex <- corridor_example(acc = TRUE)
v <- both_directions(data.frame(from = c("M03", "M08"), to = c("M20", "M25"), value = 0.05))
v$value[3:4] <- -0.05
m <- map_pairwise(v, ex$acc, ex$sites)
a <- shift_arrows(m, cells = 15)
len <- sqrt((a$xend - a$x)^2 + (a$yend - a$y)^2)
a <- a[len > 0.02 * 15 * terra::res(m)[1], ]   # too short to draw, as in plot_shift()
terra::plot(m$value, col = hcl.colors(100, "viridis"))
arrows(a$x, a$y, a$xend, a$yend, length = 0.05)

# }
```
