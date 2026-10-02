#' Transition layer from a resistance surface
#'
#' Builds the \pkg{gdistance} transition (conductance) layer used for least
#' cost paths and corridors. Conductance between neighboring cells is the
#' inverse of their mean resistance, corrected for diagonal steps. Cells at or
#' above `barrier` (for example ridgelines scored as impassable) and `NA` cells
#' cannot be crossed.
#'
#' @param resistance A single-layer `SpatRaster` (or `RasterLayer`) of
#'   resistance values, in a projected coordinate system with units of meters.
#' @param barrier Resistance at or above which a cell is a barrier. `NULL`
#'   for none.
#' @param directions 4, 8 or 16 neighbors.
#' @return A `TransitionLayer`.
#' @examples
#' ex <- corridor_example()
#' tr <- make_transition(ex$resistance)
#' @export
make_transition <- function(resistance, barrier = NULL, directions = 8) {
  r <- if (inherits(resistance, "SpatRaster")) raster::raster(resistance) else resistance
  if (!is.null(barrier)) r[r >= barrier] <- NA
  tr <- gdistance::transition(r, function(x) 1 / mean(x), directions = directions)
  gdistance::geoCorrection(tr, type = "c")
}

#' Least cost paths between sites
#'
#' @param tr A transition layer from [make_transition()].
#' @param sites An `sf` object of points or polygons (polygon centroids are
#'   used) with an ID column.
#' @param id Name of the ID column in `sites`.
#' @param pairs Optional data.frame with columns `from` and `to`. By default
#'   every unordered pair of sites.
#' @param dem Optional elevation `SpatRaster`. If given, path length is
#'   measured over the terrain surface (3D), as in Maier et al. (2022).
#' @param max_length Drop paths longer than this (meters).
#' @return An `sf` object of lines with columns `from`, `to`, `euclidean` and
#'   `length`.
#' @examples
#' ex <- corridor_example()
#' tr <- make_transition(ex$resistance)
#' paths <- least_cost_paths(tr, ex$sites, pairs = data.frame(from = "M01", to = "M12"),
#'                           dem = ex$dem)
#' paths
#' @export
least_cost_paths <- function(tr, sites, id = "site", pairs = NULL, dem = NULL, max_length = Inf) {
  xy <- site_xy(sites, id)
  if (is.null(pairs)) pairs <- all_pairs(rownames(xy))
  check_pairs(pairs, rownames(xy))
  lines <- lapply(seq_len(nrow(pairs)), function(k) {
    sl <- gdistance::shortestPath(tr, xy[pairs$from[k], ], xy[pairs$to[k], ], output = "SpatialLines")
    sf::st_geometry(sf::st_as_sf(sl))[[1]]
  })
  geom <- sf::st_sfc(lines, crs = sf::st_crs(sites))
  out <- sf::st_sf(from = pairs$from, to = pairs$to,
                   euclidean = sqrt(rowSums((xy[pairs$from, , drop = FALSE] - xy[pairs$to, , drop = FALSE])^2)),
                   geometry = geom)
  out$length <- if (is.null(dem)) as.numeric(sf::st_length(out)) else
    vapply(seq_len(nrow(out)), function(k) surface_length(geom[[k]], dem), numeric(1))
  out[out$length <= max_length, ]
}

#' Accumulated cost surfaces
#'
#' The accumulated least cost from each site to every cell of the landscape.
#' Corridors between two sites are built from the sum of their surfaces, so
#' this is computed once per site and reused for every pair.
#'
#' @inheritParams least_cost_paths
#' @param filename Optional GeoTIFF to write the surfaces to, which keeps
#'   memory use low for large landscapes.
#' @return A multi-layer `SpatRaster` with one layer per site, named by site
#'   ID. Unreachable cells are `NA`.
#' @examples
#' \donttest{
#' ex <- corridor_example()
#' tr <- make_transition(ex$resistance)
#' acc <- accumulated_cost(tr, ex$sites[1:3, ])
#' }
#' @export
accumulated_cost <- function(tr, sites, id = "site", filename = "") {
  xy <- site_xy(sites, id)
  layers <- lapply(rownames(xy), function(s) {
    a <- terra::rast(gdistance::accCost(tr, xy[s, , drop = FALSE]))
    a[is.infinite(a)] <- NA
    a
  })
  out <- terra::rast(layers)
  names(out) <- rownames(xy)
  if (nzchar(filename)) out <- terra::writeRaster(out, filename, overwrite = TRUE)
  out
}

#' Least cost corridor weights for a pair of sites
#'
#' Sums the accumulated cost surfaces of the two sites, rescales the total to
#' 0-1, keeps the cheapest `q` share of cells, and turns the kept values into
#' weights: 1 on the least cost route and 0 at the edge of the corridor.
#' Larger `q` gives a broader corridor that takes in more alternative routes.
#'
#' @param acc Accumulated cost surfaces from [accumulated_cost()].
#' @param from,to Site IDs (layer names of `acc`).
#' @param q Share of cells to keep, for example `0.001` (a narrow corridor) to
#'   `0.05` (broad). Several values give one layer each.
#' @return A `SpatRaster` of weights, `NA` outside the corridor, one layer per
#'   value of `q`.
#' @references Maier PA et al. (2022) Heredity 129:257-272.
#' @examples
#' ex <- corridor_example()
#' tr <- make_transition(ex$resistance)
#' acc <- accumulated_cost(tr, ex$sites[ex$sites$site %in% c("M01", "M12"), ])
#' w <- lcc_weights(acc, "M01", "M12", q = c(0.005, 0.05))
#' terra::plot(w)
#' @export
lcc_weights <- function(acc, from, to, q = 0.05) {
  w <- lcc_vec(acc, from, to, q)
  tmpl <- terra::rast(acc[[1]])
  out <- terra::rast(lapply(w, function(v) terra::setValues(tmpl, v)))
  names(out) <- paste0("lcc_", q)
  out
}

# Corridor weights as plain vectors (one per q), reading only the two layers
lcc_vec <- function(acc, from, to, q) {
  s <- terra::values(acc[[from]], mat = FALSE) + terra::values(acc[[to]], mat = FALSE)
  ok <- which(is.finite(s))
  v <- s[ok]
  v <- (v - min(v)) / (max(v) - min(v))
  srt <- sort(v)
  lapply(q, function(qq) {
    th <- srt[max(1L, round(length(srt) * qq))]
    w <- rep(NA_real_, length(s))
    keep <- v <= th
    # every cell on the least cost route has the same total cost; if the
    # corridor is no wider than that route, they all get weight 1
    w[ok[keep]] <- if (th > 0) 1 - v[keep] / th else 1
    w
  })
}

#' Buffers around least cost paths
#'
#' Simple corridor bandwidths: a fixed distance on either side of each least
#' cost path.
#'
#' @param paths Least cost paths from [least_cost_paths()].
#' @param width Buffer width(s) in meters.
#' @return A named list of `sf` polygon layers, one per width.
#' @examples
#' ex <- corridor_example()
#' tr <- make_transition(ex$resistance)
#' p <- least_cost_paths(tr, ex$sites, pairs = data.frame(from = "M01", to = "M12"))
#' b <- path_buffers(p, c(100, 500))
#' @export
path_buffers <- function(paths, width = c(100, 500)) {
  out <- lapply(width, function(w) sf::st_buffer(paths, w))
  names(out) <- paste0("lcp_", width, "m")
  out
}

# ---- internal helpers ---------------------------------------------------------

site_xy <- function(sites, id) {
  if (!inherits(sites, "sf")) stop("`sites` must be an sf object.", call. = FALSE)
  if (!id %in% names(sites)) stop("`sites` has no column `", id, "`.", call. = FALSE)
  g <- sf::st_geometry(sites)
  if (!all(sf::st_geometry_type(g) %in% c("POINT"))) g <- sf::st_point_on_surface(g)
  xy <- sf::st_coordinates(g)[, 1:2, drop = FALSE]
  rownames(xy) <- as.character(sites[[id]])
  if (anyDuplicated(rownames(xy))) stop("Site IDs must be unique.", call. = FALSE)
  xy
}

all_pairs <- function(ids) {
  cb <- utils::combn(ids, 2)
  data.frame(from = cb[1, ], to = cb[2, ], stringsAsFactors = FALSE)
}

check_pairs <- function(pairs, ids) {
  if (!all(c("from", "to") %in% names(pairs))) stop("`pairs` needs columns `from` and `to`.", call. = FALSE)
  miss <- setdiff(c(pairs$from, pairs$to), ids)
  if (length(miss)) stop("Unknown site IDs in `pairs`: ", paste(utils::head(miss, 5), collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

surface_length <- function(line, dem) {
  co <- sf::st_coordinates(line)[, 1:2, drop = FALSE]
  z <- terra::extract(dem, co)[, 1]
  z[is.na(z)] <- stats::approx(seq_along(z), z, seq_along(z), rule = 2)$y[is.na(z)]
  sum(sqrt(diff(co[, 1])^2 + diff(co[, 2])^2 + diff(z)^2))
}
