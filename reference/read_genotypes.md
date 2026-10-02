# Read genotypes from a file

Reads the common population genetic file formats into a
[popcounts](https://paulmaier.github.io/corridoR/reference/as_popcounts.md)
object (allele counts per population and locus), which is what
[`pairwise_fst()`](https://paulmaier.github.io/corridoR/reference/pairwise_fst.md)
and
[`directional_gst()`](https://paulmaier.github.io/corridoR/reference/directional_gst.md)
use.

## Usage

``` r
read_genotypes(
  file,
  format = c("auto", "structure", "genepop", "vcf", "genalex", "fstat", "plink"),
  pop = NULL,
  missing = c("-9"),
  onerowperind = NULL,
  label = NULL,
  popdata = NULL,
  popflag = FALSE,
  locdata = FALSE,
  phenotype = FALSE,
  extracols = 0,
  pop_names = NULL
)
```

## Arguments

- file:

  Path to the genotype file.

- format:

  One of `"auto"` (guess from the extension), `"structure"`,
  `"genepop"`, `"vcf"`, `"genalex"`, `"fstat"`, `"plink"`.

- pop:

  Optional population assignment that overrides the one in the file: a
  vector with one value per individual in file order, a two-column
  data.frame or whitespace-separated file (individual, population), or a
  function that takes individual names and returns populations (for
  example `function(x) sub("_.*", "", x)`).

- missing:

  STRUCTURE only: codes for missing alleles.

- onerowperind, label, popdata, popflag, locdata, phenotype, extracols:

  STRUCTURE only. `NULL` means detect: the row layout is detected from
  repeated labels, and the leading columns from the marker-name row when
  there is one. Otherwise the defaults are a label and a population
  column, which is the most common layout.

- pop_names:

  GENEPOP only: names for the populations, in file order. By default
  each population takes the name of its first individual, which is the
  GENEPOP convention.

## Value

A `popcounts` object.

## Supported formats

- STRUCTURE (`.str`, `.stru`, `.structure`):

  Both layouts: one row per individual (two columns per locus) or two
  rows per individual (one column per locus). Optional leading columns,
  in STRUCTURE's order: a label, the population (POPDATA), POPFLAG,
  LOCDATA, PHENOTYPE and any extra columns. An optional first row of
  marker names is detected automatically. Alleles may be any integer
  codes; missing data are `-9` by default. When the file has a
  marker-name row, the number of leading columns is worked out from it,
  otherwise set the arguments below.

- GENEPOP (`.gen`, `.genepop`):

  Title line, loci one per line or comma-separated, `Pop` separators,
  and `name , genotypes` rows. Two- or three-digit allele codes
  (detected from the genotype width), haploid or diploid; `0`, `00`,
  `000` and so on are missing.

- VCF (`.vcf`, `.vcf.gz`):

  Any ploidy, phased or unphased, biallelic or multiallelic. VCF has no
  populations, so supply `pop`.

- GenAlEx (`.csv`):

  The standard codominant layout: a first row with the number of loci,
  individuals and populations, a second row with population names, a
  header row, then sample, population and two columns per locus.

- FSTAT (`.dat`):

  The standard header (populations, loci, maximum allele, digits per
  allele), locus names, then population and genotype columns.

- PLINK (`.ped` with `.map`, or `.raw` from `--recode A`):

  The family ID is used as the population unless `pop` is given. Convert
  binary `.bed` files with `plink --recode A` or to VCF first.

Objects already in R (adegenet `genind`, `genpop`, `genlight`, or a vcfR
object) go through
[`as_popcounts()`](https://paulmaier.github.io/corridoR/reference/as_popcounts.md)
instead.

## Examples

``` r
ex <- corridor_example()
read_genotypes(ex$files$structure)
#> <popcounts> 80 populations, 1318 loci (2 alleles per locus)
#>   individuals per population: 8-8
#>   populations: M01, M02, M03, M04, M05, M06, M07, M08 ... 
read_genotypes(ex$files$genepop)
#> <popcounts> 80 populations, 1318 loci (2 alleles per locus)
#>   individuals per population: 8-8
#>   populations: M01_01, M02_01, M03_01, M04_01, M05_01, M06_01, M07_01, M08_01 ... 
read_genotypes(ex$files$vcf, pop = ex$files$popmap)
#> <popcounts> 80 populations, 1318 loci (2 alleles per locus)
#>   individuals per population: 8-8
#>   populations: M01, M02, M03, M04, M05, M06, M07, M08 ... 
```
