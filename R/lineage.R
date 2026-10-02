#' Lineage effect for pairs of sites
#'
#' Deep phylogeographic splits add genetic differentiation that has nothing to
#' do with today's landscape. Following Maier et al. (2022), the lineage
#' effect of a pair is the time to the most recent common ancestor of the two
#' sites' lineages, and zero when both belong to the same lineage. Adding it
#' as a predictor lets a model account for that history.
#'
#' @param from,to Site IDs of each pair.
#' @param lineages A data.frame with site IDs in the first column and lineage
#'   names in the second.
#' @param tmrca A data.frame with two lineage columns and the time since they
#'   split in the third (any units; only relative values matter to the
#'   models). Pairs of lineages not listed get `NA`.
#' @return A numeric vector, one value per pair.
#' @examples
#' ex <- corridor_example()
#' lin <- sf::st_drop_geometry(ex$sites)[, c("site", "lineage")]
#' lineage_cross(c("M01", "M01"), c("M02", "M30"), lin, ex$lineage_tmrca)
#' @export
lineage_cross <- function(from, to, lineages, tmrca) {
  l <- stats::setNames(as.character(lineages[[2]]), as.character(lineages[[1]]))
  miss <- setdiff(c(from, to), names(l))
  if (length(miss)) stop("No lineage for sites: ", paste(utils::head(miss, 5), collapse = ", "), call. = FALSE)
  key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "\r")
  t <- stats::setNames(as.numeric(tmrca[[3]]), key(as.character(tmrca[[1]]), as.character(tmrca[[2]])))
  a <- l[from]; b <- l[to]
  out <- unname(t[key(a, b)])
  out[a == b] <- 0
  out
}
