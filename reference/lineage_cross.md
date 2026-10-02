# Lineage effect for pairs of sites

Deep phylogeographic splits add genetic differentiation that has nothing
to do with today's landscape. Following Maier et al. (2022), the lineage
effect of a pair is the time to the most recent common ancestor of the
two sites' lineages, and zero when both belong to the same lineage.
Adding it as a predictor lets a model account for that history.

## Usage

``` r
lineage_cross(from, to, lineages, tmrca)
```

## Arguments

- from, to:

  Site IDs of each pair.

- lineages:

  A data.frame with site IDs in the first column and lineage names in
  the second.

- tmrca:

  A data.frame with two lineage columns and the time since they split in
  the third (any units; only relative values matter to the models).
  Pairs of lineages not listed get `NA`.

## Value

A numeric vector, one value per pair.

## Examples

``` r
ex <- corridor_example()
lin <- sf::st_drop_geometry(ex$sites)[, c("site", "lineage")]
lineage_cross(c("M01", "M01"), c("M02", "M30"), lin, ex$lineage_tmrca)
#> [1] 0 0
```
