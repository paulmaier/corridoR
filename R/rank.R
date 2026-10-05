#' Rank resistance surfaces by how well they explain genetic distance
#'
#' Compares competing resistance surfaces, each a hypothesis about what slows
#' movement, by how well the least cost distance through it explains a
#' pairwise genetic distance such as FST. For every surface, the cost of the
#' cheapest route between each pair of sites is related to the genetic
#' distance with a mixed model that has source and destination site as random
#' effects (Maier et al. 2022). Models are fit by maximum likelihood, so they
#' can be compared by log-likelihood and AIC. Straight-line distance is added
#' as a baseline: a surface that cannot beat it adds nothing.
#'
#' Distances are scaled to mean 0 and standard deviation 1 before fitting, so
#' the slopes are comparable across surfaces. Pairs that cannot be connected
#' through every surface (for example because a barrier isolates a site) are
#' left out of all models, so that every model is fit to the same pairs.
#'
#' @param surfaces Resistance surfaces to compare: a named list of
#'   single-layer `SpatRaster`s, or one multi-layer `SpatRaster` with a layer
#'   per hypothesis. Names label the hypotheses.
#' @param sites An `sf` object of points or polygons with an ID column.
#' @param genetic A square matrix of pairwise genetic distance (for example
#'   from [pairwise_fst()]) with site IDs as row and column names.
#' @param id Name of the ID column in `sites`.
#' @param pairs Optional data.frame with columns `from` and `to`. By default
#'   every unordered pair of sites closer than `max_dist`.
#' @param max_dist Only use pairs whose straight-line distance is below this
#'   (meters). Ignored when `pairs` is given.
#' @param barrier Resistance at or above which a cell cannot be crossed, see
#'   [make_transition()].
#' @param straight Add straight-line distance as a baseline hypothesis.
#' @return A data.frame with one row per hypothesis, sorted from best to
#'   worst: `hypothesis`, `logLik`, `AIC`, `delta_AIC` (difference from the
#'   best model), `R2m` (marginal R2, the share of variance explained by
#'   distance alone; Nakagawa & Schielzeth 2013) and `slope` (the effect of
#'   one standard deviation of distance). Attribute `distances` holds the
#'   pairs and their distance under each hypothesis.
#' @references
#'   Maier PA, Vandergast AG, Ostoja SM, Aguilar A, Bohonak AJ (2022) Landscape
#'   genetics of a sub-alpine toad: climate change predicted to induce upward
#'   range shifts via asymmetrical migration corridors. Heredity 129:257-272.
#'
#'   Nakagawa S, Schielzeth H (2013) A general and simple method for obtaining
#'   R2 from generalized linear mixed-effects models. Methods in Ecology and
#'   Evolution 4:133-142.
#' @examples
#' \donttest{
#' ex <- corridor_example()
#' fst <- pairwise_fst(read_genotypes(ex$files$structure))
#' surfaces <- list(slope = ex$resistance_slope, cover = ex$resistance_cover,
#'                  slope_80_cover_20 = ex$resistance)
#' rank_resistance(surfaces, ex$sites, fst, max_dist = 15000)
#' }
#' @export
rank_resistance <- function(surfaces, sites, genetic, id = "site", pairs = NULL, max_dist = Inf,
                            barrier = 1e6, straight = TRUE) {
  need("lme4")
  if (inherits(surfaces, "SpatRaster")) {
    surfaces <- stats::setNames(lapply(seq_len(terra::nlyr(surfaces)), function(k) surfaces[[k]]),
                                names(surfaces))
  }
  if (!is.list(surfaces) || is.null(names(surfaces)) || any(names(surfaces) == "") ||
      anyDuplicated(names(surfaces))) {
    stop("`surfaces` must be a named list of SpatRasters (or a multi-layer SpatRaster) with unique names.",
         call. = FALSE)
  }
  xy <- site_xy(sites, id)
  if (is.null(rownames(genetic)) || !all(rownames(xy) %in% rownames(genetic))) {
    stop("`genetic` needs row and column names matching the site IDs.", call. = FALSE)
  }
  euclid <- function(p) sqrt(rowSums((xy[p$from, , drop = FALSE] - xy[p$to, , drop = FALSE])^2))
  if (is.null(pairs)) {
    pairs <- all_pairs(rownames(xy))
    pairs <- pairs[euclid(pairs) < max_dist, , drop = FALSE]
  }
  check_pairs(pairs, rownames(xy))

  dist <- vapply(surfaces, function(r) {
    cd <- as.matrix(gdistance::costDistance(make_transition(r, barrier = barrier), xy))
    dimnames(cd) <- list(rownames(xy), rownames(xy))
    cd[cbind(pairs$from, pairs$to)]
  }, numeric(nrow(pairs)))
  dist <- matrix(dist, nrow = nrow(pairs), dimnames = list(NULL, names(surfaces)))
  if (straight) dist <- cbind(dist, straight_line = euclid(pairs))

  y <- genetic[cbind(pairs$from, pairs$to)]
  ok <- is.finite(y) & apply(is.finite(dist), 1, all)
  if (any(!ok)) {
    warning(sum(!ok), " of ", length(ok), " pairs could not be connected through every surface ",
            "(or lack a genetic distance) and were left out.", call. = FALSE)
  }
  d <- data.frame(from = factor(pairs$from[ok]), to = factor(pairs$to[ok]), y = y[ok])

  res <- lapply(colnames(dist), function(h) {
    d$x <- as.numeric(scale(dist[ok, h]))
    fit <- suppressMessages(lme4::lmer(y ~ x + (1 | from) + (1 | to), data = d, REML = FALSE))
    v_fixed <- stats::var(as.vector(lme4::getME(fit, "X") %*% lme4::fixef(fit)))
    v_random <- sum(as.numeric(lme4::VarCorr(fit)))
    data.frame(hypothesis = h, logLik = as.numeric(stats::logLik(fit)), AIC = stats::AIC(fit),
               R2m = v_fixed / (v_fixed + v_random + stats::sigma(fit)^2),
               slope = unname(lme4::fixef(fit)["x"]))
  })
  out <- do.call(rbind, res)
  out <- out[order(out$AIC), , drop = FALSE]
  out$delta_AIC <- out$AIC - out$AIC[1]
  out <- out[, c("hypothesis", "logLik", "AIC", "delta_AIC", "R2m", "slope")]
  rownames(out) <- NULL
  attr(out, "distances") <- data.frame(pairs[ok, c("from", "to")], genetic = y[ok], dist[ok, , drop = FALSE],
                                       row.names = NULL)
  out
}
