---
title: 'corridoR: Migration corridors and genetic range shifts from landscape genetic data'
tags:
  - R
  - landscape genetics
  - least cost corridors
  - gene flow
  - climate change
  - range shifts
authors:
  - name: Paul A. Maier
    orcid: 0000-0003-0851-8827
    affiliation: "1, 2"
affiliations:
  - name: FamilyTreeDNA, Gene by Gene, Houston, Texas, USA
    index: 1
  - name: Department of Biology, San Diego State University, San Diego, California, USA
    index: 2
date: 1 October 2026
bibliography: paper.bib
---

<!-- FIRST DRAFT for the author to rewrite. Facts, structure and references
     are in place; the wording is placeholder. JOSS asks for 250-1000 words. -->

# Summary

Climate change is expected to move species upslope and poleward, but for most
species we cannot watch individuals disperse. Genetic data record where
migrants have come from and gone to. `corridoR` turns pairwise genetic
differentiation and directional migration among sampled sites into maps of
how connectivity, and the net direction of migration, are projected to change
under future climate. It implements the workflow introduced by @maier2022,
which forecast upslope range shifts in the Yosemite toad (*Anaxyrus canorus*).

The package covers the steps from genotype files to maps. It reads the common
population genetic formats, computes pairwise FST and the directional
differentiation of @sundqvist2016, builds least cost corridors between sites,
summarizes raw environmental rasters inside each corridor weighted by how
likely a route through each cell is, chooses a corridor width for each group
of environmental features by cross-validation, fits Cubist models
[@quinlan1992; @kuhn2023] of FST and net migration, and maps projected changes
back onto the landscape as magnitude and net-direction surfaces.

# Statement of need

Landscape genetic studies usually encode each hypothesis as a single
resistance surface and compare surfaces by how well resistance distance
explains genetic distance [@mcrae2006; @peterman2018]. This works for one or
two variables, but translating raw values into resistance is subjective, and
the number of combinations grows quickly as variables are added. Transect
methods extract raw values along straight lines between sites
[@vanstrien2012], which keeps the variables raw but ignores the routes animals
actually take.

`corridoR` separates the two problems. A simple, well-supported resistance
model defines the likely routes, and least cost corridors around them become
the units over which raw environmental features are summarized. Each feature
group can choose its own bandwidth, so broad-scale influences such as climate
and local ones such as soil are each measured at a fitting scale. Because the
features stay raw, future climate layers can be swapped in directly to
forecast change. Mapping pairwise forecasts back onto their corridors, with
direction vectors for asymmetric migration, gives a spatial picture of where
connectivity is projected to fall and which way net movement will turn.

The package also restores tools that are hard to find in current R. The
directional differentiation of @sundqvist2016 was available through
`diveRsity::divMigrate`, which has since been archived on CRAN; `corridoR`
reimplements it in plain R and reads STRUCTURE, GENEPOP, VCF, GenAlEx, FSTAT
and PLINK files as well as `adegenet` [@jombart2008] and `vcfR` [@knaus2017]
objects.

# State of the field

`gdistance` [@vanetten2017] computes least cost paths and accumulated costs,
`ResistanceGA` [@peterman2018] optimizes resistance surfaces, and Circuitscape
[@mcrae2006] models current flow. None of these summarize raw covariates
within weighted corridors, select corridor widths per feature group, or map
pairwise forecasts of directional migration. `corridoR` builds on `gdistance`,
`terra` [@hijmans2024] and `sf` [@pebesma2018] and is meant to sit alongside
these tools rather than replace them.

# Functionality

- `read_genotypes()` and `as_popcounts()`: genotype input; `pairwise_fst()`
  and `directional_gst()`: genetic distances and net migration.
- `make_transition()`, `least_cost_paths()`, `accumulated_cost()`,
  `lcc_weights()` and `path_buffers()`: routes and corridors.
- `corridor_extract()` and `site_contrast()`: between-site and at-site
  features.
- `select_bandwidth()`, `fit_connectivity()` and `predict()`: models and
  forecasts.
- `map_pairwise()`, `shift_arrows()` and `plot_shift()`: maps of change and net direction.
- `lineage_cross()`: the phylogeographic lineage effect.

A simulated mountain range ships with the package: 80 meadows in four
lineages on a gradual western slope, a crest and a steep eastern escarpment,
with warming that dries low meadows the most. The tutorial runs the full
workflow on it and recovers the upslope shift built into the simulation, with
the direction of the largest projected changes correct for 98% of pairs.

# AI usage

Claude (Anthropic) helped with fast implementation and hypothesis-testing,
while the design, the validation, and every scientific decision stay with me.

# Acknowledgements

<!-- Funding and thanks; the 2022 study was funded by the U.S. Geological
     Survey Natural Resource Preservation Program. -->

# References
