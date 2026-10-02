# Case weights for pairwise data

Sites that appear in many pairs would otherwise dominate a model of
pairwise data. Each pair is weighted by one minus the mean relative
frequency of its two sites, then the weights are shifted so the largest
is 1 and normalized to sum to 1, as in Maier et al. (2022).

## Usage

``` r
pair_weights(from, to)
```

## Arguments

- from, to:

  Site IDs of each pair.

## Value

Numeric weights, one per pair.

## Examples

``` r
pair_weights(c("a", "a", "a", "b"), c("b", "c", "d", "c"))
#> [1] 0.2272727 0.2272727 0.2727273 0.2727273
```
