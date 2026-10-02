#' Simulated example landscape
#'
#' Loads the example shipped with the package: a simulated 60 x 40 km
#' mountain range modeled loosely on the central Sierra Nevada, with 80
#' snowmelt meadows in four lineages. A long, gradual western slope rises to a
#' crest near the eastern edge; east of it a steep escarpment lies in a rain
#' shadow; two steep-walled canyons cut the western slope. Migrants leave
#' water-stressed meadows for wetter neighbors. Under warming (+4.2 C)
#' snowpack shrinks by a share that falls with elevation, from about 90% in
#' the foothills to under 20% at the crest, so the push uphill grows most at
#' low elevation. The script that builds the example is in the package source
#' (`data-raw/make_example.R`).
#'
#' @return A list with
#' \describe{
#' \item{`dem`, `resistance`}{Elevation (m) and movement resistance
#'   `SpatRaster`s; ridgelines have resistance 1e6.}
#' \item{`env`, `env_future`}{Environmental layers now and under warming:
#'   snowpack (mm), runoff (mm), a meadow moisture index (log scale), summer
#'   temperature (C), slope (degrees) and forest cover (0-1).}
#' \item{`sites`}{The 80 meadows as `sf` points with `site`, `elevation` and
#'   `lineage`.}
#' \item{`lineage_tmrca`}{Time since each pair of lineages split (thousand
#'   years), for [lineage_cross()].}
#' \item{`resistance_slope`, `resistance_cover`}{The two components blended
#'   80:20 into `resistance`, for comparing path hypotheses.}
#' \item{`acc`}{Accumulated cost surfaces for every meadow, from
#'   [accumulated_cost()]. Only with `acc = TRUE`; they take a few seconds to
#'   compute and are cached for the rest of the session.}
#' \item{`files`}{Paths to the genotypes (about 1,300 SNPs, 8 individuals per meadow)
#'   in STRUCTURE (`structure`), GENEPOP (`genepop`) and VCF (`vcf`) format,
#'   the VCF population map (`popmap`), and a small GENEPOP sample of six
#'   meadows and 100 SNPs (`small`) for quick tests.}
#' }
#' @param acc Also compute the accumulated cost surfaces.
#' @examples
#' ex <- corridor_example()
#' terra::plot(ex$dem)
#' plot(sf::st_geometry(ex$sites), add = TRUE, pch = 19)
#' @export
corridor_example <- function(acc = FALSE) {
  f <- function(x) system.file("extdata", x, package = "corridoR", mustWork = TRUE)
  land <- terra::rast(f("landscape.tif"))
  out <- list(dem = land[["elevation"]], resistance = land[["resistance"]],
       resistance_slope = land[["resistance_slope"]], resistance_cover = land[["resistance_cover"]],
       lineage_tmrca = utils::read.csv(f("lineage_tmrca.csv")),
       env = terra::rast(f("env_present.tif")), env_future = terra::rast(f("env_future.tif")),
       sites = sf::st_read(f("sites.gpkg"), quiet = TRUE),
       files = list(structure = f("example.str"), genepop = f("example.gen"),
                    vcf = f("example.vcf"), popmap = f("popmap.txt"),
                    small = f("example_small.gen")))
  if (acc) {
    if (is.null(.cache$acc)) {
      tr <- make_transition(out$resistance, barrier = 1e6)
      file <- file.path(tempdir(), "corridoR_example_acc.tif")
      .cache$acc <- accumulated_cost(tr, out$sites, filename = file)
    }
    out$acc <- .cache$acc
  }
  out
}

.cache <- new.env(parent = emptyenv())
