# Pairwise FST between populations

Pairwise FST between populations

## Usage

``` r
pairwise_fst(x, method = c("hudson", "nei"))
```

## Arguments

- x:

  A
  [popcounts](https://paulmaier.github.io/corridoR/reference/as_popcounts.md)
  object, or anything
  [`as_popcounts()`](https://paulmaier.github.io/corridoR/reference/as_popcounts.md)
  accepts.

- method:

  `"hudson"` (default) is Hudson's FST as a ratio of averages across
  loci (Bhatia et al. 2013), which is robust to unequal sample sizes.
  `"nei"` is Nei's GST from pooled allele frequencies.

## Value

A symmetric matrix of pairwise FST with population names as dimnames.
Any precomputed matrix (for example from hierfstat or Stacks) works just
as well in the rest of the package.

## References

Bhatia G, Patterson N, Sankararaman S, Price AL (2013) Estimating and
interpreting FST: the impact of rare variants. Genome Research
23:1514-1521.

## Examples

``` r
ex <- corridor_example()
fst <- pairwise_fst(read_genotypes(ex$files$structure))
round(fst[1:4, 1:4], 3)
#>       M01   M02   M03   M04
#> M01 0.000 0.163 0.154 0.240
#> M02 0.163 0.000 0.166 0.241
#> M03 0.154 0.166 0.000 0.235
#> M04 0.240 0.241 0.235 0.000
```
