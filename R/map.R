#' Map pairwise values onto the landscape
#'
#' Spreads one value per pair of sites (a predicted FST, a projected change in
#' dM, and so on) over the pair's least cost corridor, scaled by the corridor
#' weight so that less likely routes count for less, and averages overlapping
#' corridors cell by cell. With `direction = TRUE`, each directed pair also
#' contributes a vector from source to destination, reversed when its value is
#' negative, and the vectors are summed in every cell to give the net
#' direction of change, as in Fig. 1D of Maier et al. (2022).
#'
#' @param values A data.frame with columns `from`, `to` and `value`. For
#'   directional maps include both directions of each pair.
#' @param acc Accumulated cost surfaces from [accumulated_cost()].
#' @param sites An `sf` object of sites (used for the direction of each pair).
#' @param id Name of the ID column in `sites`.
#' @param q Corridor threshold passed to [lcc_weights()].
#' @param direction Also compute the net direction of change.
#' @param vectors How each pair contributes to the direction. `"value"` (the
#'   default) adds a unit vector from source to destination scaled by the
#'   pair's value, so pairs with larger changes pull harder. `"sign"` adds the
#'   source-to-destination displacement scaled only by the sign of the value,
#'   as in Maier et al. (2022).
#' @param absolute Map the absolute value (the magnitude of change), as for
#'   change maps. Set to `FALSE` to map signed values.
#' @param progress Show a progress bar.
#' @return A `SpatRaster` with layer `value` (mean weighted value per cell) and,
#'   with `direction = TRUE`, `vx` and `vy` (mean direction vector, map units)
#'   and `bearing` (degrees clockwise from north).
#' @examples
#' \donttest{
#' ex <- corridor_example(acc = TRUE)
#' # a made-up change for a few pairs: more net movement toward the east
#' v <- both_directions(data.frame(from = c("M03", "M08", "M14"), to = c("M20", "M25", "M28")))
#' east <- sf::st_coordinates(ex$sites)[match(v$to, ex$sites$site), 1] >
#'   sf::st_coordinates(ex$sites)[match(v$from, ex$sites$site), 1]
#' v$value <- ifelse(east, 0.05, -0.05)
#' m <- map_pairwise(v, ex$acc, ex$sites)
#' terra::plot(m$value, col = hcl.colors(100, "viridis"))
#' }
#' @export
map_pairwise <- function(values, acc, sites, id = "site", q = 0.05, direction = TRUE,
                         vectors = c("value", "sign"), absolute = TRUE, progress = interactive()) {
  vectors <- match.arg(vectors)
  if (!all(c("from", "to", "value") %in% names(values))) {
    stop("`values` needs columns `from`, `to` and `value`.", call. = FALSE)
  }
  xy <- site_xy(sites, id)
  check_pairs(values, rownames(xy))
  n <- terra::ncell(acc)
  sum_v <- cnt <- vx <- vy <- numeric(n)
  key <- paste(pmin(values$from, values$to), pmax(values$from, values$to))
  groups <- split(seq_len(nrow(values)), key)
  pb <- if (progress) utils::txtProgressBar(0, length(groups), style = 3) else NULL
  for (k in seq_along(groups)) {
    rows <- groups[[k]]
    w <- lcc_vec(acc, values$from[rows[1]], values$to[rows[1]], q)[[1]]
    on <- which(!is.na(w) & w > 0)
    w <- w[on]
    for (r in rows) {
      v <- values$value[r]
      if (is.na(v)) next
      sum_v[on] <- sum_v[on] + w * (if (absolute) abs(v) else v)
      cnt[on] <- cnt[on] + 1
      if (direction && v != 0) {
        d <- xy[values$to[r], ] - xy[values$from[r], ]
        # "value": unit direction scaled by the change; "sign": the pair's
        # displacement scaled by the sign of the change (Maier et al. 2022)
        d <- if (vectors == "value") d / sqrt(sum(d^2)) * v else d * sign(v)
        vx[on] <- vx[on] + w * d[1]
        vy[on] <- vy[on] + w * d[2]
      }
    }
    if (!is.null(pb)) utils::setTxtProgressBar(pb, k)
  }
  if (!is.null(pb)) close(pb)
  avg <- function(s) { o <- s / cnt; o[cnt == 0] <- NA; o }
  tmpl <- terra::rast(acc[[1]])
  out <- terra::setValues(tmpl, avg(sum_v))
  names(out) <- "value"
  if (direction) {
    ax <- avg(vx); ay <- avg(vy)
    out <- c(out, terra::setValues(tmpl, ax), terra::setValues(tmpl, ay),
             terra::setValues(tmpl, (atan2(ax, ay) * 180 / pi) %% 360))
    names(out) <- c("value", "vx", "vy", "bearing")
  }
  out
}

#' Arrows for a directional map
#'
#' Summarizes a map from [map_pairwise()] on a coarser grid, for drawing
#' arrows. Direction follows the summed vectors. By default length does too,
#' so arrows are long where corridors agree on a direction and shrink where
#' opposing flows cancel; `length_by = "magnitude"` scales them by the mean
#' size of change instead.
#'
#' @param map A directional map from [map_pairwise()].
#' @param cells Grid size, in cells of `map`.
#' @param length Length of the longest arrow, in map units. Defaults to
#'   about 80% of a grid cell.
#' @param length_by `"net"` (length of the summed vector) or `"magnitude"`
#'   (mean size of change in the grid cell).
#' @param scale `"linear"` lengths, or `"sqrt"` to compress the range so a
#'   few very large changes do not shrink every other arrow to a dot.
#' @param min_quantile Drop arrows in the weakest cells (below this quantile
#'   of magnitude).
#' @param min_length Shortest arrow, as a share of `length`. With the default
#'   of 0, arrow length is proportional to magnitude, so arrows shrink to
#'   nothing where nothing changes.
#' @return A data.frame with `x`, `y`, `xend`, `yend` and `magnitude`, ready
#'   for `ggplot2::geom_segment()` or [graphics::arrows()].
#' @examples
#' \donttest{
#' ex <- corridor_example(acc = TRUE)
#' v <- both_directions(data.frame(from = c("M03", "M08"), to = c("M20", "M25"), value = 0.05))
#' v$value[3:4] <- -0.05
#' m <- map_pairwise(v, ex$acc, ex$sites)
#' a <- shift_arrows(m, cells = 15)
#' len <- sqrt((a$xend - a$x)^2 + (a$yend - a$y)^2)
#' a <- a[len > 0.02 * 15 * terra::res(m)[1], ]   # too short to draw, as in plot_shift()
#' terra::plot(m$value, col = hcl.colors(100, "viridis"))
#' arrows(a$x, a$y, a$xend, a$yend, length = 0.05)
#' }
#' @export
shift_arrows <- function(map, cells = 25, length = NULL, length_by = c("net", "magnitude"),
                         scale = c("linear", "sqrt"), min_quantile = 0, min_length = 0) {
  length_by <- match.arg(length_by)
  scale <- match.arg(scale)
  mag <- terra::aggregate(map$value, cells, mean, na.rm = TRUE)
  sx <- terra::aggregate(map$vx, cells, sum, na.rm = TRUE)
  sy <- terra::aggregate(map$vy, cells, sum, na.rm = TRUE)
  d <- as.data.frame(c(mag, sx, sy), xy = TRUE, na.rm = TRUE)
  names(d) <- c("x", "y", "magnitude", "vx", "vy")
  d <- d[d$magnitude >= stats::quantile(d$magnitude, min_quantile) & (d$vx != 0 | d$vy != 0), ]
  if (is.null(length)) length <- 0.8 * cells * terra::res(map)[1]
  h <- sqrt(d$vx^2 + d$vy^2)
  size <- if (length_by == "net") h else d$magnitude
  if (scale == "sqrt") size <- sqrt(size)
  len <- length * (min_length + (1 - min_length) * size / max(size))
  data.frame(x = d$x, y = d$y, xend = d$x + d$vx / h * len, yend = d$y + d$vy / h * len,
             magnitude = d$magnitude)
}

#' Plot a map of projected change
#'
#' Draws a map from [map_pairwise()] over shaded relief. Change is colored on
#' a three-tone scale: cool blue below `threshold`, a narrow white band at it,
#' and red above it, so the areas under most pressure stand out. With a
#' directional map, arrows from [shift_arrows()] show the net direction.
#'
#' @param map A map from [map_pairwise()].
#' @param dem Optional elevation `SpatRaster` for the relief.
#' @param sites Optional `sf` points to draw on top.
#' @param col_sites Fill colors for `sites` (recycled, or one per site).
#' @param threshold Value drawn in white. Defaults to the median of the map.
#' @param smooth Width (in cells) of a moving-average window applied to the
#'   displayed surface only; 1 for none.
#' @param alpha Opacity of the change layer.
#' @param arrows Draw arrows (needs a directional map).
#' @param cells,scale Passed to [shift_arrows()].
#' @param main,legend_title Title and legend title.
#' @return The map, invisibly.
#' @examples
#' \donttest{
#' ex <- corridor_example(acc = TRUE)
#' v <- both_directions(data.frame(from = c("M03", "M08"), to = c("M40", "M45"), value = 0.05))
#' v$value[3:4] <- -0.05
#' plot_shift(map_pairwise(v, ex$acc, ex$sites), dem = ex$dem, sites = ex$sites)
#' }
#' @export
plot_shift <- function(map, dem = NULL, sites = NULL, col_sites = "white", threshold = NULL,
                       smooth = 5, alpha = 0.5, arrows = "vx" %in% names(map), cells = 14,
                       scale = c("linear", "sqrt"), main = "", legend_title = "|change|") {
  v <- map$value
  if (smooth > 1) v <- terra::focal(v, w = smooth, fun = mean, na.rm = TRUE, na.policy = "omit")
  if (is.null(threshold)) threshold <- stats::median(terra::values(v, mat = FALSE), na.rm = TRUE)
  mx <- terra::global(v, "max", na.rm = TRUE)[[1]]
  br <- c(seq(0, threshold * 0.92, length.out = 20), threshold, seq(threshold * 1.08, mx, length.out = 20))
  cols <- c(grDevices::colorRampPalette(c("#5b8fc7", "#bcd6ee"))(19), "#f7f7f7", "#f7f7f7",
            grDevices::colorRampPalette(c("#f4a582", "#b2182b", "#67001f"))(19))
  if (!is.null(dem)) {
    hill <- terra::shade(terra::terrain(dem * 3, "slope", unit = "radians"),
                         terra::terrain(dem, "aspect", unit = "radians"), 35, 315)
    terra::plot(hill, col = grDevices::grey(0:100 / 100), legend = FALSE, axes = FALSE,
                mar = c(0.5, 0.5, 2, 6), main = main)
    tint <- grDevices::colorRampPalette(c("#b5a06e", "#cdbf8f", "#dfe0c4", "#f2f2ee", "#ffffff"))(60)
    terra::plot(dem, col = tint, alpha = 0.35, add = TRUE, legend = FALSE)
    terra::plot(v, col = cols, breaks = br, alpha = alpha, add = TRUE, type = "continuous",
                plg = list(title = legend_title))
  } else {
    terra::plot(v, col = cols, breaks = br, alpha = alpha, type = "continuous", axes = FALSE,
                mar = c(0.5, 0.5, 2, 6), main = main, plg = list(title = legend_title))
  }
  if (arrows) {
    a <- shift_arrows(map, cells = cells, scale = match.arg(scale))
    long <- sqrt((a$xend - a$x)^2 + (a$yend - a$y)^2) > 0.02 * cells * terra::res(map)[1]
    graphics::points(a$x[!long], a$y[!long], pch = 16, cex = 0.25)
    a <- a[long, ]
    graphics::arrows(a$x, a$y, a$xend, a$yend, length = 0.05, lwd = 1.5)
  }
  if (!is.null(sites)) {
    graphics::points(sf::st_coordinates(sites)[, 1:2, drop = FALSE], pch = 21, bg = col_sites,
                     col = "white", cex = 1.4, lwd = 1)
  }
  invisible(map)
}
