# Changelog

## corridoR (development version)

- [`fit_connectivity()`](https://paulmaier.github.io/corridoR/reference/fit_connectivity.md)
  and
  [`select_bandwidth()`](https://paulmaier.github.io/corridoR/reference/select_bandwidth.md)
  now keep both directions of a pair (A to B and B to A) in the same
  cross-validation fold. Before, a pair could be predicted from its own
  reverse, which inflated cross-validated R2 and could favor overfit
  settings.

## corridoR 0.1.0

- First release: genotype readers (STRUCTURE, GENEPOP, VCF, GenAlEx,
  FSTAT, PLINK, adegenet and vcfR objects), pairwise FST and directional
  GST, least cost paths and corridors, weighted environmental
  extraction, per-group bandwidth selection, Cubist models of
  connectivity, and maps of projected change in magnitude and net
  direction. Includes a simulated example landscape.
