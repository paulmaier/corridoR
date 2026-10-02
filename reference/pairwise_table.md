# Pairwise genetic table

Turns matrices of FST and dM into a long table with one row per directed
pair, the layout used by the modeling functions.

## Usage

``` r
pairwise_table(fst, dM = NULL)
```

## Arguments

- fst:

  Symmetric matrix of pairwise differentiation.

- dM:

  Optional matrix of net migration, `dM[i, j]` from `i` to `j`.

## Value

A data.frame with columns `from`, `to`, `fst` and (if given) `dM`.

## Examples

``` r
ex <- corridor_example()
pc <- read_genotypes(ex$files$structure)
head(pairwise_table(pairwise_fst(pc), directional_gst(pc)$dM))
#>   from  to        fst           dM
#> 1  M02 M01 0.16311080  0.038685611
#> 2  M03 M01 0.15357827  0.004226984
#> 3  M04 M01 0.23965549  0.011150260
#> 4  M05 M01 0.08890528 -0.088469572
#> 5  M06 M01 0.22114611  0.013136436
#> 6  M07 M01 0.13162280 -0.027356787
```
