# Transition layer from a resistance surface

Builds the gdistance transition (conductance) layer used for least cost
paths and corridors. Conductance between neighboring cells is the
inverse of their mean resistance, corrected for diagonal steps. Cells at
or above `barrier` (for example ridgelines scored as impassable) and
`NA` cells cannot be crossed.

## Usage

``` r
make_transition(resistance, barrier = NULL, directions = 8)
```

## Arguments

- resistance:

  A single-layer `SpatRaster` (or `RasterLayer`) of resistance values,
  in a projected coordinate system with units of meters.

- barrier:

  Resistance at or above which a cell is a barrier. `NULL` for none.

- directions:

  4, 8 or 16 neighbors.

## Value

A `TransitionLayer`.

## Examples

``` r
# \donttest{
ex <- corridor_example()
tr <- make_transition(ex$resistance)
# }
```
