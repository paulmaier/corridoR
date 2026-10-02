#' Choose a corridor bandwidth for each group of features
#'
#' For each feature group, fits a random forest of the response (usually
#' FST) on that group's features, once per candidate bandwidth, and keeps the
#' bandwidth with the lowest cross-validated RMSE. This lets broad-scale
#' influences such as climate pick a wide corridor while local ones such as
#' soil pick a narrow one.
#'
#' @param features A named list of data.frames, one per bandwidth, each with
#'   columns `from`, `to` and the features, as returned by
#'   [corridor_extract()]. All must list the same pairs in the same order.
#' @param response Numeric response for each row (each pair).
#' @param groups A named list mapping group names to feature (column) names.
#' @param weights Optional case weights, see [pair_weights()].
#' @param folds Number of cross-validation folds.
#' @param num.trees Trees per forest.
#' @param min.node.size Values to try for the minimum node size. `mtry` is
#'   tuned over every value from 1 to the number of features in the group.
#' @param seed Random seed for the folds and forests.
#' @return A data.frame with one row per group and bandwidth (group,
#'   bandwidth, mtry, min.node.size, RMSE, Rsquared), sorted by group and
#'   RMSE, with attribute `best`: a named vector of the chosen bandwidth per
#'   group.
#' @examples
#' # Two candidate bandwidths; climate measured in the wide corridor explains
#' # the response better
#' set.seed(1)
#' n <- 150
#' wide <- data.frame(from = "a", to = "b", snow = rnorm(n), slope = rnorm(n))
#' narrow <- data.frame(from = "a", to = "b", snow = wide$snow + rnorm(n), slope = rnorm(n))
#' fst <- 0.1 - 0.03 * wide$snow + rnorm(n, sd = 0.01)
#' sel <- select_bandwidth(list(narrow = narrow, wide = wide), fst,
#'                         groups = list(climate = "snow", terrain = "slope"),
#'                         num.trees = 100, folds = 3)
#' attr(sel, "best")
#' @export
select_bandwidth <- function(features, response, groups, weights = NULL, folds = 5,
                             num.trees = 1000, min.node.size = c(1, 5, 10), seed = 12345) {
  need("ranger")
  set.seed(seed)
  fold_id <- sample(rep(seq_len(folds), length.out = length(response)))
  res <- list()
  for (g in names(groups)) for (bw in names(features)) {
    x <- features[[bw]][, groups[[g]], drop = FALSE]
    x <- x[, vapply(x, function(v) length(unique(v)) > 1, logical(1)), drop = FALSE]
    if (!ncol(x)) next
    grid <- expand.grid(mtry = seq_len(ncol(x)), min.node.size = min.node.size)
    sc <- t(vapply(seq_len(nrow(grid)), function(k) {
      pred <- numeric(length(response))
      for (f in seq_len(folds)) {
        tr <- fold_id != f
        fit <- ranger::ranger(x = x[tr, , drop = FALSE], y = response[tr], num.trees = num.trees,
                              mtry = grid$mtry[k], min.node.size = grid$min.node.size[k],
                              case.weights = if (is.null(weights)) NULL else weights[tr],
                              seed = seed, num.threads = 1)
        pred[!tr] <- stats::predict(fit, x[!tr, , drop = FALSE])$predictions
      }
      c(rmse = sqrt(mean((pred - response)^2)), r2 = stats::cor(pred, response)^2)
    }, numeric(2)))
    b <- which.min(sc[, "rmse"])
    res[[length(res) + 1]] <- data.frame(group = g, bandwidth = bw, mtry = grid$mtry[b],
                                         min.node.size = grid$min.node.size[b],
                                         RMSE = sc[b, "rmse"], Rsquared = sc[b, "r2"])
  }
  out <- do.call(rbind, res)
  out <- out[order(out$group, out$RMSE), ]
  rownames(out) <- NULL
  best <- vapply(split(out, out$group), function(d) d$bandwidth[which.min(d$RMSE)], character(1))
  attr(out, "best") <- best
  out
}

#' Case weights for pairwise data
#'
#' Sites that appear in many pairs would otherwise dominate a model of
#' pairwise data. Each pair is weighted by one minus the mean relative
#' frequency of its two sites, then the weights are shifted so the largest is
#' 1 and normalized to sum to 1, as in Maier et al. (2022).
#'
#' @param from,to Site IDs of each pair.
#' @return Numeric weights, one per pair.
#' @examples
#' pair_weights(c("a", "a", "a", "b"), c("b", "c", "d", "c"))
#' @export
pair_weights <- function(from, to) {
  freq <- table(c(from, to)); freq <- freq / max(freq)
  w <- 1 - (as.numeric(freq[as.character(from)]) + as.numeric(freq[as.character(to)])) / 2
  w <- w + (1 - max(w))
  w / sum(w)
}

#' Fit a Cubist model of connectivity
#'
#' The modeling step of Maier et al. (2022). Features are first reduced by a
#' PCA within each group (groups of more than two features, centered and
#' scaled), which removes most redundancy among related layers. A Cubist
#' model (rule-based trees with a linear model in each leaf, which
#' extrapolates better than random forests) is tuned by cross-validation over
#' committees and neighbors. Its variable importance then decides which
#' feature to keep whenever features from different groups are collinear
#' (VIF above `vif_threshold`), and the model is refit on what remains.
#'
#' @param data A data.frame of features, one row per directed pair.
#' @param response Numeric response (FST, dM, or another pairwise statistic).
#' @param groups Named list mapping group names to columns of `data`. Columns
#'   not in any group (for example path length or a lineage term) are used
#'   as they are and always kept.
#' @param committees,neighbors Values to try.
#' @param folds Number of cross-validation folds.
#' @param vif_threshold Collinearity cutoff. `Inf` skips the VIF step.
#' @param seed Random seed for the folds.
#' @return An object of class `corridor_model` with the final Cubist fit,
#'   the PCA for each group, the features kept, the cross-validation results
#'   and predictions, and variable importance. Use [stats::predict()] for new or
#'   future data.
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 300
#' d <- data.frame(snow = rnorm(n), runoff = rnorm(n), temp = rnorm(n), dist = runif(n))
#' d$runoff <- d$snow + rnorm(n, sd = 0.3)
#' y <- 0.05 + 0.02 * d$snow + 0.03 * d$dist + rnorm(n, sd = 0.01)
#' m <- fit_connectivity(d, y, groups = list(climate = c("snow", "runoff", "temp")),
#'                       committees = c(1, 10), neighbors = c(0, 5), folds = 5)
#' m
#' future <- transform(d, snow = snow - 1)
#' summary(predict(m, future) - predict(m, d))
#' }
#' @export
fit_connectivity <- function(data, response, groups, committees = c(1, 10, 50, 75, 100),
                             neighbors = c(0, 1, 5, 7, 9), folds = 10, vif_threshold = 10,
                             seed = 12345) {
  need("Cubist")
  pcs <- lapply(names(groups), function(g) group_pca(data[, groups[[g]], drop = FALSE], g))
  names(pcs) <- names(groups)
  fixed <- setdiff(names(data), c(unlist(groups), "from", "to"))
  X <- cbind(data[, fixed, drop = FALSE], do.call(cbind, unname(lapply(pcs, `[[`, "x"))))
  X <- X[, vapply(X, function(v) length(unique(v)) > 1, logical(1)), drop = FALSE]

  full <- tune_cubist(X, response, committees, neighbors, folds, seed)
  vars <- names(X)
  if (is.finite(vif_threshold)) {
    imp <- cubist_importance(full$fit)
    cand <- setdiff(names(X), fixed)
    keep <- vif_select(X[, cand, drop = FALSE], imp[cand], vif_threshold)
    vars <- c(intersect(fixed, names(X)), keep)
  }
  final <- tune_cubist(X[, vars, drop = FALSE], response, committees, neighbors, folds, seed)
  imp <- cubist_importance(final$fit)
  structure(list(fit = final$fit, committees = final$best$committees, neighbors = final$best$neighbors,
                 vars = vars, groups = groups, fixed = fixed,
                 pca = lapply(pcs, `[[`, "pca"), cv = final$results, best = final$best,
                 cv_pred = final$cv_pred, observed = response,
                 importance = sort(100 * imp / max(imp), decreasing = TRUE)),
            class = "corridor_model")
}

#' @export
print.corridor_model <- function(x, ...) {
  cat(sprintf("<corridor_model> Cubist, %d committees, %d neighbors\n", x$committees, x$neighbors))
  cat(sprintf("  %d features after PCA and VIF filtering; CV RMSE %.4g, R2 %.3f\n",
              length(x$vars), x$best$RMSE, x$best$Rsquared))
  top <- utils::head(x$importance, 5)
  cat("  most important:", paste(sprintf("%s (%.0f)", names(top), top), collapse = ", "), "\n")
  invisible(x)
}

#' Predict from a connectivity model
#'
#' Projects new feature values (for example future climate) onto the PCA
#' loadings fitted to present-day data, so the components keep their meaning,
#' then predicts with the Cubist model.
#'
#' @param object A `corridor_model` from [fit_connectivity()].
#' @param newdata A data.frame with the same feature columns as the data the
#'   model was fit to.
#' @param ... Not used.
#' @return Numeric predictions.
#' @export
predict.corridor_model <- function(object, newdata, ...) {
  parts <- list(newdata[, intersect(object$fixed, names(newdata)), drop = FALSE])
  for (g in names(object$groups)) {
    x <- newdata[, object$groups[[g]], drop = FALSE]
    p <- object$pca[[g]]
    parts[[g]] <- if (is.null(p)) x else {
      s <- as.data.frame(stats::predict(p, x[, rownames(p$rotation), drop = FALSE]))
      stats::setNames(s, paste0(g, "_", names(s)))
    }
  }
  X <- do.call(cbind, unname(parts))
  stats::predict(object$fit, X[, object$vars, drop = FALSE], neighbors = object$neighbors)
}

# ---- internal helpers ---------------------------------------------------------

need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("This function needs the '", pkg, "' package: install.packages(\"", pkg, "\")", call. = FALSE)
  }
}

group_pca <- function(x, prefix) {
  x <- x[, vapply(x, function(v) length(unique(v)) > 1, logical(1)), drop = FALSE]
  if (ncol(x) <= 2) return(list(x = x, pca = NULL))
  pca <- stats::prcomp(x, scale. = TRUE)
  s <- as.data.frame(pca$x)
  names(s) <- paste0(prefix, "_", names(s))
  list(x = s, pca = pca)
}

tune_cubist <- function(X, y, committees, neighbors, folds, seed) {
  set.seed(seed)
  fold_id <- sample(rep(seq_len(folds), length.out = length(y)))
  pred <- array(NA_real_, c(length(y), length(committees), length(neighbors)))
  for (f in seq_len(folds)) {
    tr <- fold_id != f
    for (ci in seq_along(committees)) {
      fit <- Cubist::cubist(X[tr, , drop = FALSE], y[tr], committees = committees[ci])
      for (ni in seq_along(neighbors)) {
        pred[!tr, ci, ni] <- stats::predict(fit, X[!tr, , drop = FALSE], neighbors = neighbors[ni])
      }
    }
  }
  grid <- expand.grid(committees = committees, neighbors = neighbors)
  grid$RMSE <- grid$Rsquared <- NA_real_
  for (k in seq_len(nrow(grid))) {
    p <- pred[, match(grid$committees[k], committees), match(grid$neighbors[k], neighbors)]
    grid$RMSE[k] <- sqrt(mean((p - y)^2)); grid$Rsquared[k] <- stats::cor(p, y)^2
  }
  b <- which.min(grid$RMSE)
  fit <- Cubist::cubist(X, y, committees = grid$committees[b])
  list(fit = fit, results = grid, best = grid[b, ],
       cv_pred = pred[, match(grid$committees[b], committees), match(grid$neighbors[b], neighbors)])
}

cubist_importance <- function(fit) {
  u <- fit$usage
  imp <- stats::setNames((u$Conditions + u$Model) / 2, u$Variable)
  out <- stats::setNames(numeric(length(fit$vars$all)), fit$vars$all)
  out[names(imp)] <- imp
  out
}

# Drop features until every VIF is below the threshold. Features are added in
# order of importance; one is kept only if it does not push any VIF over the
# threshold.
vif_select <- function(X, importance, threshold) {
  ord <- names(sort(importance[names(X)], decreasing = TRUE))
  keep <- character(0)
  for (v in ord) {
    trial <- c(keep, v)
    if (length(trial) < 2 || max(vif(X[, trial, drop = FALSE])) <= threshold) keep <- trial
  }
  keep
}

vif <- function(X) {
  X <- as.matrix(X)
  r <- stats::cor(X)
  out <- tryCatch(diag(solve(r)), error = function(e) rep(Inf, ncol(X)))
  stats::setNames(out, colnames(X))
}
