# corridoR (development version)

* The example data no longer raise warnings with older GDAL and PROJ
  libraries: the sites are stored as GeoPackage 1.2, and the example
  coordinate system names its datum (WGS84) explicitly. Coordinates and
  values are unchanged.

* New `rank_resistance()` compares candidate resistance surfaces by how well
  least cost distance through each explains genetic distance (mixed models
  with source and destination as random effects), with straight-line distance
  as a baseline. The tutorial uses it in step 2.

* `fit_connectivity()` and `select_bandwidth()` now keep both directions of a
  pair (A to B and B to A) in the same cross-validation fold. Before, a pair
  could be predicted from its own reverse, which inflated cross-validated R2
  and could favor overfit settings.

# corridoR 0.1.0

* First release: genotype readers (STRUCTURE, GENEPOP, VCF, GenAlEx, FSTAT,
  PLINK, adegenet and vcfR objects), pairwise FST and directional GST, least
  cost paths and corridors, weighted environmental extraction, per-group
  bandwidth selection, Cubist models of connectivity, and maps of projected
  change in magnitude and net direction. Includes a simulated example
  landscape.
