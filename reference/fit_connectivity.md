# Fit a Cubist model of connectivity

The modeling step of Maier et al. (2022). Features are first reduced by
a PCA within each group (groups of more than two features, centered and
scaled), which removes most redundancy among related layers. A Cubist
model (rule-based trees with a linear model in each leaf, which
extrapolates better than random forests) is tuned by cross-validation
over committees and neighbors. Its variable importance then decides
which feature to keep whenever features from different groups are
collinear (VIF above `vif_threshold`), and the model is refit on what
remains.

## Usage

``` r
fit_connectivity(
  data,
  response,
  groups,
  committees = c(1, 10, 50, 75, 100),
  neighbors = c(0, 1, 5, 7, 9),
  folds = 10,
  vif_threshold = 10,
  seed = 12345
)
```

## Arguments

- data:

  A data.frame of features, one row per directed pair.

- response:

  Numeric response (FST, dM, or another pairwise statistic).

- groups:

  Named list mapping group names to columns of `data`. Columns not in
  any group (for example path length or a lineage term) are used as they
  are and always kept.

- committees, neighbors:

  Values to try.

- folds:

  Number of cross-validation folds. When `data` has `from` and `to`
  columns, both directions of a pair always fall in the same fold, so a
  pair is never used to predict its own reverse.

- vif_threshold:

  Collinearity cutoff. `Inf` skips the VIF step.

- seed:

  Random seed for the folds.

## Value

An object of class `corridor_model` with the final Cubist fit, the PCA
for each group, the features kept, the cross-validation results and
predictions, and variable importance. Use
[`stats::predict()`](https://rdrr.io/r/stats/predict.html) for new or
future data.

## Examples

``` r
# \donttest{
set.seed(1)
n <- 300
d <- data.frame(snow = rnorm(n), runoff = rnorm(n), temp = rnorm(n), dist = runif(n))
d$runoff <- d$snow + rnorm(n, sd = 0.3)
y <- 0.05 + 0.02 * d$snow + 0.03 * d$dist + rnorm(n, sd = 0.01)
m <- fit_connectivity(d, y, groups = list(climate = c("snow", "runoff", "temp")),
                      committees = c(1, 10), neighbors = c(0, 5), folds = 5)
m
#> <corridor_model> Cubist, 10 committees, 0 neighbors
#>   4 features after PCA and VIF filtering; CV RMSE 0.01019, R2 0.833
#>   most important: dist (100), climate_PC1 (100), climate_PC3 (100), climate_PC2 (0) 
future <- transform(d, snow = snow - 1)
summary(predict(m, future) - predict(m, d))
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -0.02086 -0.02086 -0.02086 -0.02086 -0.02086 -0.02086 
# }
```
