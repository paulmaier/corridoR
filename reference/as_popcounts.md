# Population allele counts

Every genetic statistic in corridoR works from allele counts per
population and locus. A `popcounts` object holds them: a list with one
integer matrix per locus (populations in rows, alleles in columns), plus
the population names and the number of individuals sampled in each.

## Usage

``` r
as_popcounts(x, ...)

# S3 method for class 'popcounts'
as_popcounts(x, ...)

# S3 method for class 'genind'
as_popcounts(x, pop = NULL, ...)

# S3 method for class 'genpop'
as_popcounts(x, ...)

# S3 method for class 'genlight'
as_popcounts(x, pop = NULL, ...)

# S3 method for class 'vcfR'
as_popcounts(x, pop = NULL, ...)
```

## Arguments

- x:

  An object to convert.

- ...:

  Not used.

- pop:

  Population of each individual, for objects that do not carry it (a
  `vcfR` object, or a `genlight` without `pop`). See
  [`read_genotypes()`](https://paulmaier.github.io/corridoR/reference/read_genotypes.md)
  for the accepted forms.

## Value

An object of class `popcounts`.

## Details

`as_popcounts()` converts objects from other packages: `genind`,
`genpop` and `genlight` from adegenet and `vcfR` from vcfR. To read a
file directly, use
[`read_genotypes()`](https://paulmaier.github.io/corridoR/reference/read_genotypes.md).

## Examples

``` r
ex <- corridor_example()
pc <- read_genotypes(ex$files$structure)
pc
#> <popcounts> 80 populations, 1318 loci (2 alleles per locus)
#>   individuals per population: 8-8
#>   populations: M01, M02, M03, M04, M05, M06, M07, M08 ... 
```
