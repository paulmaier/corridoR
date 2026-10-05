# corridoR

**Migration corridors and genetic range shifts from landscape genetic
data.**

corridoR models genetic connectivity along migration corridors and
forecasts how it shifts under climate change. It implements the workflow
of [Maier et al. (2022,
*Heredity*)](https://doi.org/10.1038/s41437-022-00561-x), which
predicted upslope range shifts in the Yosemite toad:

1.  **Genetic distances.** Read genotypes from STRUCTURE, GENEPOP, VCF,
    GenAlEx, FSTAT or PLINK files (or `adegenet` and `vcfR` objects),
    then compute pairwise FST and directional differentiation (dM, the
    net direction of migration).
2.  **Least cost corridors.** Widen least cost paths into corridors
    whose cells are weighted by how likely a route through them is.
3.  **Raw environmental features.** Summarize any raster inside each
    corridor, plus the contrast between source and destination sites.
4.  **A bandwidth for each feature group.** Let climate act over a broad
    corridor and soil over a narrow one, chosen by cross-validation.
5.  **Cubist models** of FST and dM, which extrapolate to future climate
    better than random forests.
6.  **Maps** of projected change, with arrows showing the net direction
    of migration.

![The simulated example: 80 meadows in four lineages on a range that
rises from west to east](reference/figures/study_area.png)

The simulated example: 80 meadows in four lineages on a range that rises
from west to east

![Projected shift in net migration: low meadows are pushed upslope, the
crest is a refuge](reference/figures/shift_map.png)

Projected shift in net migration: low meadows are pushed upslope, the
crest is a refuge

## Installation

``` r

# install.packages("remotes")
remotes::install_github("paulmaier/corridoR")
```

## A quick look

``` r

library(corridoR)
ex <- corridor_example(acc = TRUE)          # a simulated mountain range, 80 meadows

pc  <- read_genotypes(ex$files$structure)   # or GENEPOP, VCF, ...
fst <- pairwise_fst(pc)
dM  <- directional_gst(pc)$dM

# corridor weights between two meadows, and climate summarized inside them
w <- lcc_weights(ex$acc, "M20", "M60", q = 0.05)
corridor_extract(ex$env, data.frame(from = "M20", to = "M60"), acc = ex$acc, q = 0.05)
```

The [full
tutorial](https://paulmaier.github.io/corridoR/articles/tutorial.html)
runs all six steps on the simulated range, from genotype files to the
maps above, and the [getting-started
vignette](https://paulmaier.github.io/corridoR/articles/corridoR.html)
is a short version. The same analysis on real data is in the [Yosemite
toad
tutorial](https://www.paulmaierresearch.com/software/yosemite-toad-corridors/).

## Genotype formats

| Format | Read with | Notes |
|----|----|----|
| STRUCTURE (`.str`) | [`read_genotypes()`](https://paulmaier.github.io/corridoR/reference/read_genotypes.md) | one or two rows per individual; optional label, POPDATA, POPFLAG, LOCDATA, PHENOTYPE and extra columns; marker-name row detected |
| GENEPOP (`.gen`) | [`read_genotypes()`](https://paulmaier.github.io/corridoR/reference/read_genotypes.md) | 2- or 3-digit alleles, haploid or diploid |
| VCF (`.vcf`, `.vcf.gz`) | `read_genotypes(pop = )` | any ploidy, phased or unphased, multiallelic |
| GenAlEx (`.csv`) | [`read_genotypes()`](https://paulmaier.github.io/corridoR/reference/read_genotypes.md) | codominant layout |
| FSTAT (`.dat`) | [`read_genotypes()`](https://paulmaier.github.io/corridoR/reference/read_genotypes.md) |  |
| PLINK (`.ped`/`.map`, `.raw`) | [`read_genotypes()`](https://paulmaier.github.io/corridoR/reference/read_genotypes.md) | convert `.bed` with `plink --recode A` |
| `genind`, `genpop`, `genlight`, `vcfR` | [`as_popcounts()`](https://paulmaier.github.io/corridoR/reference/as_popcounts.md) | from adegenet or vcfR |

Precomputed matrices of FST or dM from any other software work too.

## Citation

If you use corridoR, **please cite both papers**:

> Maier PA, Vandergast AG, Ostoja SM, Aguilar A, Bohonak AJ (2022)
> Landscape genetics of a sub-alpine toad: climate change predicted to
> induce upward range shifts via asymmetrical migration corridors.
> *Heredity* 129:257-272. <https://doi.org/10.1038/s41437-022-00561-x>

> Maier PA (2027) corridoR: an R package for forecasting genetic range
> shifts along migration corridors. Manuscript in preparation.

The first introduced the method; the second describes the package.
`citation("corridoR")` prints both.

## License

MIT
