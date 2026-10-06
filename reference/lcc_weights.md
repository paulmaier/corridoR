# Least cost corridor weights for a pair of sites

Sums the accumulated cost surfaces of the two sites, rescales the total
to 0-1, keeps the cheapest `q` share of cells, and turns the kept values
into weights: 1 on the least cost route and 0 at the edge of the
corridor. Larger `q` gives a broader corridor that takes in more
alternative routes.

## Usage

``` r
lcc_weights(acc, from, to, q = 0.05)
```

## Arguments

- acc:

  Accumulated cost surfaces from
  [`accumulated_cost()`](https://paulmaier.github.io/corridoR/reference/accumulated_cost.md).

- from, to:

  Site IDs (layer names of `acc`).

- q:

  Share of cells to keep, for example `0.001` (a narrow corridor) to
  `0.05` (broad). Several values give one layer each.

## Value

A `SpatRaster` of weights, `NA` outside the corridor, one layer per
value of `q`.

## References

Maier PA et al. (2022) Heredity 129:257-272.

## Examples

``` r
# \donttest{
ex <- corridor_example()
tr <- make_transition(ex$resistance)
acc <- accumulated_cost(tr, ex$sites[ex$sites$site %in% c("M01", "M12"), ])
w <- lcc_weights(acc, "M01", "M12", q = c(0.005, 0.05))
terra::plot(w, col = hcl.colors(100, "viridis"))

# }
```
