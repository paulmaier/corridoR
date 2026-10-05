# Rank resistance surfaces by how well they explain genetic distance

Compares competing resistance surfaces, each a hypothesis about what
slows movement, by how well the least cost distance through it explains
a pairwise genetic distance such as FST. For every surface, the cost of
the cheapest route between each pair of sites is related to the genetic
distance with a mixed model that has source and destination site as
random effects (Maier et al. 2022). Models are fit by maximum
likelihood, so they can be compared by log-likelihood and AIC.
Straight-line distance is added as a baseline: a surface that cannot
beat it adds nothing.

## Usage

``` r
rank_resistance(
  surfaces,
  sites,
  genetic,
  id = "site",
  pairs = NULL,
  max_dist = Inf,
  barrier = 1e+06,
  straight = TRUE
)
```

## Arguments

- surfaces:

  Resistance surfaces to compare: a named list of single-layer
  `SpatRaster`s, or one multi-layer `SpatRaster` with a layer per
  hypothesis. Names label the hypotheses.

- sites:

  An `sf` object of points or polygons with an ID column.

- genetic:

  A square matrix of pairwise genetic distance (for example from
  [`pairwise_fst()`](https://paulmaier.github.io/corridoR/reference/pairwise_fst.md))
  with site IDs as row and column names.

- id:

  Name of the ID column in `sites`.

- pairs:

  Optional data.frame with columns `from` and `to`. By default every
  unordered pair of sites closer than `max_dist`.

- max_dist:

  Only use pairs whose straight-line distance is below this (meters).
  Ignored when `pairs` is given.

- barrier:

  Resistance at or above which a cell cannot be crossed, see
  [`make_transition()`](https://paulmaier.github.io/corridoR/reference/make_transition.md).

- straight:

  Add straight-line distance as a baseline hypothesis.

## Value

A data.frame with one row per hypothesis, sorted from best to worst:
`hypothesis`, `logLik`, `AIC`, `delta_AIC` (difference from the best
model), `R2m` (marginal R2, the share of variance explained by distance
alone; Nakagawa & Schielzeth 2013) and `slope` (the effect of one
standard deviation of distance). Attribute `distances` holds the pairs
and their distance under each hypothesis.

## Details

Distances are scaled to mean 0 and standard deviation 1 before fitting,
so the slopes are comparable across surfaces. Pairs that cannot be
connected through every surface (for example because a barrier isolates
a site) are left out of all models, so that every model is fit to the
same pairs.

## References

Maier PA, Vandergast AG, Ostoja SM, Aguilar A, Bohonak AJ (2022)
Landscape genetics of a sub-alpine toad: climate change predicted to
induce upward range shifts via asymmetrical migration corridors.
Heredity 129:257-272.

Nakagawa S, Schielzeth H (2013) A general and simple method for
obtaining R2 from generalized linear mixed-effects models. Methods in
Ecology and Evolution 4:133-142.

## Examples

``` r
# \donttest{
ex <- corridor_example()
fst <- pairwise_fst(read_genotypes(ex$files$structure))
surfaces <- list(slope = ex$resistance_slope, cover = ex$resistance_cover,
                 slope_80_cover_20 = ex$resistance)
rank_resistance(surfaces, ex$sites, fst, max_dist = 15000)
#>          hypothesis   logLik       AIC delta_AIC        R2m       slope
#> 1             slope 2288.989 -4567.978   0.00000 0.03622447 0.009864794
#> 2 slope_80_cover_20 2275.221 -4540.442  27.53579 0.05508055 0.011909024
#> 3             cover 2121.487 -4232.974 335.00414 0.10355020 0.016687517
#> 4     straight_line 2115.328 -4220.656 347.32251 0.01264270 0.006066638
# }
```
