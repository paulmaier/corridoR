# Use each pair in both directions

Corridor features are the same whichever site is the source, so they are
measured once per unordered pair. Models of directional statistics such
as dM need a row for each direction; this repeats every row with `from`
and `to` swapped.

## Usage

``` r
both_directions(x)
```

## Arguments

- x:

  A data.frame with columns `from` and `to`.

## Value

`x` with the reversed pairs appended.

## Examples

``` r
both_directions(data.frame(from = "a", to = "b", snow = 3))
#>   from to snow
#> 1    a  b    3
#> 2    b  a    3
```
