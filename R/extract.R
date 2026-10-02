#' Summarize environmental layers within corridors
#'
#' For every pair of sites, summarizes each environmental layer inside the
#' pair's corridor. With least cost corridors (`acc` and `q`), each cell
#' counts in proportion to its corridor weight, so likely routes matter more
#' than marginal ones. With path buffers, cells are weighted by the share of
#' the cell inside the buffer.
#'
#' @param env A `SpatRaster` of environmental layers. Layers are resampled to
#'   the corridor grid as needed (bilinear for continuous layers).
#' @param pairs Data.frame with columns `from` and `to`.
#' @param acc,q Accumulated cost surfaces from [accumulated_cost()] and the
#'   corridor threshold for [lcc_weights()]. Use these for least cost
#'   corridors.
#' @param buffers Alternatively, an `sf` layer of buffered paths from
#'   [path_buffers()] with columns `from` and `to`.
#' @param stat How to summarize each layer: `"mean"` (weighted mean, the
#'   default) or `"sum"` (weighted sum, for counts such as stream or trail
#'   crossings). A named vector sets the statistic per layer.
#' @param progress Show a progress bar.
#' @return A data.frame with `from`, `to` and one column per layer of `env`.
#' @examples
#' ex <- corridor_example()
#' tr <- make_transition(ex$resistance)
#' pr <- data.frame(from = c("M01", "M02"), to = c("M12", "M05"))
#' acc <- accumulated_cost(tr, ex$sites[ex$sites$site %in% unlist(pr), ])
#' corridor_extract(ex$env, pr, acc = acc, q = 0.05)
#' @export
corridor_extract <- function(env, pairs, acc = NULL, q = 0.05, buffers = NULL,
                             stat = "mean", progress = interactive()) {
  check_env(env)
  stat <- rep_stat(stat, names(env))
  if (is.null(acc) == is.null(buffers)) stop("Give either `acc` (least cost corridors) or `buffers`.", call. = FALSE)
  pb <- if (progress) utils::txtProgressBar(0, nrow(pairs), style = 3) else NULL
  if (!is.null(acc)) {
    env_c <- align_env(env, acc)
    ev <- terra::values(env_c)                         # cells x layers, read once
    is_sum <- vapply(stat, identical, logical(1), "sum")
    rows <- lapply(seq_len(nrow(pairs)), function(k) {
      w <- lcc_vec(acc, pairs$from[k], pairs$to[k], q[1])[[1]]
      on <- which(!is.na(w))
      if (!is.null(pb)) utils::setTxtProgressBar(pb, k)
      weighted_cols(ev[on, , drop = FALSE], w[on], is_sum)
    })
  } else {
    key <- paste(buffers$from, buffers$to)
    rows <- lapply(seq_len(nrow(pairs)), function(k) {
      i <- match(paste(pairs$from[k], pairs$to[k]), key)
      if (is.na(i)) i <- match(paste(pairs$to[k], pairs$from[k]), key)
      if (is.na(i)) stop("No buffer for pair ", pairs$from[k], "-", pairs$to[k], call. = FALSE)
      x <- terra::extract(env, terra::vect(buffers[i, ]), exact = TRUE)
      if (!is.null(pb)) utils::setTxtProgressBar(pb, k)
      vapply(names(env), function(v) {
        ok <- !is.na(x[[v]])
        s <- sum(x[[v]][ok] * x$fraction[ok])
        if (stat[[v]] == "sum") s else s / sum(x$fraction[ok])
      }, numeric(1))
    })
  }
  if (!is.null(pb)) close(pb)
  out <- as.data.frame(do.call(rbind, rows))
  cbind(data.frame(from = pairs$from, to = pairs$to, stringsAsFactors = FALSE), out)
}

#' Environmental contrasts between sites
#'
#' At-site features: the value of each environmental layer at the source site
#' minus its value at the destination site. Sites can be points (the cell
#' value is used) or polygons (the area-weighted mean is used). Site-level
#' attributes such as meadow area or network degree can be supplied as a
#' data.frame instead of, or in addition to, rasters.
#'
#' @param sites An `sf` object of sites.
#' @param pairs Data.frame with columns `from` and `to` (directed pairs).
#' @param env Optional `SpatRaster` of environmental layers.
#' @param attributes Optional data.frame with an ID column and one column per
#'   site attribute.
#' @param id Name of the ID column in `sites` (and `attributes`).
#' @param suffix Appended to the feature names.
#' @return A data.frame with `from`, `to` and one column per feature, each
#'   equal to the source value minus the destination value.
#' @examples
#' ex <- corridor_example()
#' pr <- data.frame(from = c("M01", "M12"), to = c("M12", "M01"))
#' site_contrast(ex$sites, pr, env = ex$env[["snowpack"]])
#' @export
site_contrast <- function(sites, pairs, env = NULL, attributes = NULL, id = "site", suffix = ".at") {
  ids <- as.character(sites[[id]])
  check_pairs(pairs, ids)
  vals <- data.frame(row.names = ids)
  if (!is.null(env)) {
    check_env(env)
    g <- sf::st_geometry(sites)
    if (all(sf::st_geometry_type(g) == "POINT")) {
      v <- terra::extract(env, terra::vect(sites), ID = FALSE)
    } else {
      x <- terra::extract(env, terra::vect(sites), exact = TRUE)
      v <- do.call(rbind, lapply(split(x, x$ID), function(d) {
        vapply(names(env), function(n) stats::weighted.mean(d[[n]], d$fraction, na.rm = TRUE), numeric(1))
      }))
    }
    vals <- cbind(vals, as.data.frame(v))
  }
  if (!is.null(attributes)) {
    a <- attributes[match(ids, as.character(attributes[[id]])), setdiff(names(attributes), id), drop = FALSE]
    vals <- cbind(vals, a)
  }
  d <- as.matrix(vals[pairs$from, , drop = FALSE]) - as.matrix(vals[pairs$to, , drop = FALSE])
  colnames(d) <- paste0(colnames(vals), suffix)
  cbind(data.frame(from = pairs$from, to = pairs$to, stringsAsFactors = FALSE), as.data.frame(d, row.names = NULL))
}

# ---- internal helpers ---------------------------------------------------------

check_env <- function(env) {
  if (!inherits(env, "SpatRaster")) stop("`env` must be a terra SpatRaster.", call. = FALSE)
  if (anyDuplicated(names(env))) stop("Layer names in `env` must be unique.", call. = FALSE)
}

rep_stat <- function(stat, layers) {
  if (length(stat) == 1 && is.null(names(stat))) stat <- stats::setNames(rep(stat, length(layers)), layers)
  stat <- stat[layers]
  stat[is.na(stat)] <- "mean"
  if (!all(stat %in% c("mean", "sum"))) stop("`stat` must be \"mean\" or \"sum\".", call. = FALSE)
  as.list(stats::setNames(stat, layers))
}

align_env <- function(env, template) {
  if (terra::compareGeom(env, template[[1]], stopOnError = FALSE)) return(env)
  terra::resample(env, template[[1]], method = "bilinear")
}

# Weighted sums or means of each column; cells missing in a layer are skipped
weighted_cols <- function(m, w, is_sum) {
  ok <- !is.na(m)
  m[!ok] <- 0
  s <- colSums(m * w)
  out <- ifelse(is_sum, s, s / colSums(ok * w))
  stats::setNames(out, colnames(m))
}

#' Use each pair in both directions
#'
#' Corridor features are the same whichever site is the source, so they are
#' measured once per unordered pair. Models of directional statistics such as
#' dM need a row for each direction; this repeats every row with `from` and
#' `to` swapped.
#'
#' @param x A data.frame with columns `from` and `to`.
#' @return `x` with the reversed pairs appended.
#' @examples
#' both_directions(data.frame(from = "a", to = "b", snow = 3))
#' @export
both_directions <- function(x) {
  r <- x
  r$from <- x$to; r$to <- x$from
  out <- rbind(x, r)
  rownames(out) <- NULL
  out
}
