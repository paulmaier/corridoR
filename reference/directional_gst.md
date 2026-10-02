# Directional genetic differentiation and net migration (dM)

Computes the directional GST of Sundqvist et al. (2016), as in
`diveRsity::divMigrate(stat = "gst")`: for each pair of populations, a
hypothetical shared gene pool is built from the normalized geometric
mean of their allele frequencies, and each population's differentiation
from that pool is converted to relative migration. Following Maier et
al. (2022), the difference between emigration and immigration gives net
migration, dM.

## Usage

``` r
directional_gst(x)
```

## Arguments

- x:

  A
  [popcounts](https://paulmaier.github.io/corridoR/reference/as_popcounts.md)
  object, or anything
  [`as_popcounts()`](https://paulmaier.github.io/corridoR/reference/as_popcounts.md)
  accepts.

## Value

A list with

- `relative`:

  Relative migration, scaled so the largest value is 1. `relative[i, j]`
  is migration from population `i` into population `j`.

- `dM`:

  `relative - t(relative)`: positive when `i` sends more migrants to `j`
  than it receives from it.

- `gst`:

  The directional GST matrix behind `relative`.

## References

Sundqvist L, Keenan K, Zackrisson M, Prodohl P, Kleinhans D (2016)
Directional genetic differentiation and relative migration. Ecology and
Evolution 6:3461-3475.

Maier PA, Vandergast AG, Ostoja SM, Aguilar A, Bohonak AJ (2022)
Landscape genetics of a sub-alpine toad: climate change predicted to
induce upward range shifts via asymmetrical migration corridors.
Heredity 129:257-272.

## Examples

``` r
ex <- corridor_example()
dg <- directional_gst(read_genotypes(ex$files$structure))
round(dg$dM[1:4, 1:4], 3)
#>       M01    M02    M03    M04
#> M01 0.000 -0.039 -0.004 -0.011
#> M02 0.039  0.000  0.031  0.001
#> M03 0.004 -0.031  0.000 -0.032
#> M04 0.011 -0.001  0.032  0.000
```
