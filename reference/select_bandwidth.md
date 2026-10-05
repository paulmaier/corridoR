# Choose a corridor bandwidth for each group of features

For each feature group, fits a random forest of the response (usually
FST) on that group's features, once per candidate bandwidth, and keeps
the bandwidth with the lowest cross-validated RMSE. This lets
broad-scale influences such as climate pick a wide corridor while local
ones such as soil pick a narrow one.

## Usage

``` r
select_bandwidth(
  features,
  response,
  groups,
  weights = NULL,
  folds = 5,
  num.trees = 1000,
  min.node.size = c(1, 5, 10),
  seed = 12345
)
```

## Arguments

- features:

  A named list of data.frames, one per bandwidth, each with columns
  `from`, `to` and the features, as returned by
  [`corridor_extract()`](https://paulmaier.github.io/corridoR/reference/corridor_extract.md).
  All must list the same pairs in the same order.

- response:

  Numeric response for each row (each pair).

- groups:

  A named list mapping group names to feature (column) names.

- weights:

  Optional case weights, see
  [`pair_weights()`](https://paulmaier.github.io/corridoR/reference/pair_weights.md).

- folds:

  Number of cross-validation folds. Both directions of a pair (rows A to
  B and B to A) always fall in the same fold.

- num.trees:

  Trees per forest.

- min.node.size:

  Values to try for the minimum node size. `mtry` is tuned over every
  value from 1 to the number of features in the group.

- seed:

  Random seed for the folds and forests.

## Value

A data.frame with one row per group and bandwidth (group, bandwidth,
mtry, min.node.size, RMSE, Rsquared), sorted by group and RMSE, with
attribute `best`: a named vector of the chosen bandwidth per group.

## Examples

``` r
# Two candidate bandwidths; climate measured in the wide corridor explains
# the response better
set.seed(1)
n <- 150
ids <- data.frame(from = paste0("s", 1:n), to = paste0("t", 1:n))
wide <- data.frame(ids, snow = rnorm(n), slope = rnorm(n))
narrow <- data.frame(ids, snow = wide$snow + rnorm(n), slope = rnorm(n))
fst <- 0.1 - 0.03 * wide$snow + rnorm(n, sd = 0.01)
sel <- select_bandwidth(list(narrow = narrow, wide = wide), fst,
                        groups = list(climate = "snow", terrain = "slope"),
                        num.trees = 100, folds = 3)
attr(sel, "best")
#> climate terrain 
#>  "wide"  "wide" 
```
