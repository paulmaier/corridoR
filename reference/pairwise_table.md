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
fst <- matrix(c(0, 0.1, 0.2, 0.1, 0, 0.15, 0.2, 0.15, 0), 3,
              dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
pairwise_table(fst)
#>   from to  fst
#> 1    b  a 0.10
#> 2    c  a 0.20
#> 3    a  b 0.10
#> 4    c  b 0.15
#> 5    a  c 0.20
#> 6    b  c 0.15
```
