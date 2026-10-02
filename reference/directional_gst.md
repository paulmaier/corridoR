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
dg <- directional_gst(read_genotypes(ex$files$small))
round(dg$dM[1:4, 1:4], 3)
#>        M03_01 M15_01 M28_01 M41_01
#> M03_01  0.000 -0.021  0.099 -0.101
#> M15_01  0.021  0.000 -0.258 -0.043
#> M28_01 -0.099  0.258  0.000  0.290
#> M41_01  0.101  0.043 -0.290  0.000
```
