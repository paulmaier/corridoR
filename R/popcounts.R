#' Population allele counts
#'
#' Every genetic statistic in corridoR works from allele counts per
#' population and locus. A `popcounts` object holds them: a list with one
#' integer matrix per locus (populations in rows, alleles in columns), plus the
#' population names and the number of individuals sampled in each.
#'
#' `as_popcounts()` converts objects from other packages: `genind`, `genpop`
#' and `genlight` from \pkg{adegenet} and `vcfR` from \pkg{vcfR}. To read a file
#' directly, use [read_genotypes()].
#'
#' @param x An object to convert.
#' @param pop Population of each individual, for objects that do not carry it
#'   (a `vcfR` object, or a `genlight` without `pop`). See [read_genotypes()]
#'   for the accepted forms.
#' @param ... Not used.
#' @return An object of class `popcounts`.
#' @examples
#' ex <- corridor_example()
#' pc <- read_genotypes(ex$files$small)
#' pc
#' @export
as_popcounts <- function(x, ...) UseMethod("as_popcounts")

#' @rdname as_popcounts
#' @export
as_popcounts.popcounts <- function(x, ...) x

#' @rdname as_popcounts
#' @export
as_popcounts.genind <- function(x, pop = NULL, ...) {
  tab <- x@tab
  p <- resolve_pop(if (is.null(pop)) x@pop else pop, rownames(tab))
  loc <- as.character(x@loc.fac)
  counts <- lapply(split(seq_len(ncol(tab)), factor(loc, levels = unique(loc))), function(j) {
    m <- tab[, j, drop = FALSE]
    m[is.na(m)] <- 0L
    out <- rowsum(m, p, reorder = FALSE)
    colnames(out) <- sub("^[^.]*\\.", "", colnames(tab)[j])
    out[levels(p), , drop = FALSE]
  })
  new_popcounts(counts, levels(p), table(p)[levels(p)])
}

#' @rdname as_popcounts
#' @export
as_popcounts.genpop <- function(x, ...) {
  tab <- x@tab
  loc <- as.character(x@loc.fac)
  counts <- lapply(split(seq_len(ncol(tab)), factor(loc, levels = unique(loc))), function(j) {
    m <- tab[, j, drop = FALSE]
    m[is.na(m)] <- 0L
    m
  })
  new_popcounts(counts, rownames(tab), rep(NA_integer_, nrow(tab)))
}

#' @rdname as_popcounts
#' @export
as_popcounts.genlight <- function(x, pop = NULL, ...) {
  g <- as.matrix(x)
  p <- resolve_pop(if (is.null(pop)) x@pop else pop, rownames(g))
  ploidy <- if (length(x@ploidy)) max(x@ploidy) else 2L
  dosage_counts(g, p, ploidy)
}

#' @rdname as_popcounts
#' @export
as_popcounts.vcfR <- function(x, pop = NULL, ...) {
  gt <- x@gt[, -1, drop = FALSE]
  fmt <- x@gt[, 1]
  gt_field <- vapply(strsplit(fmt, ":", fixed = TRUE), function(f) match("GT", f), integer(1))
  calls <- gt
  for (i in seq_len(nrow(gt))) calls[i, ] <- vapply(strsplit(gt[i, ], ":", fixed = TRUE), `[`, "", gt_field[i])
  ids <- x@fix[, "ID"]
  ids[is.na(ids) | ids == "."] <- paste(x@fix[, "CHROM"], x@fix[, "POS"], sep = "_")[is.na(ids) | ids == "."]
  vcf_gt_counts(calls, ids, resolve_pop(pop, colnames(gt)))
}

new_popcounts <- function(counts, pops, n_ind) {
  counts <- lapply(counts, function(m) {
    storage.mode(m) <- "integer"
    m <- m[, colSums(m) > 0, drop = FALSE]
    rownames(m) <- pops
    m
  })
  keep <- vapply(counts, ncol, integer(1)) > 0
  structure(list(counts = counts[keep], pops = as.character(pops),
                 n_ind = stats::setNames(as.integer(n_ind), pops)),
            class = "popcounts")
}

#' @export
print.popcounts <- function(x, ...) {
  na <- vapply(x$counts, ncol, integer(1))
  cat(sprintf("<popcounts> %d populations, %d loci (%s alleles per locus)\n",
              length(x$pops), length(x$counts),
              if (length(unique(na)) == 1) na[1] else paste(range(na), collapse = "-")))
  if (!all(is.na(x$n_ind))) {
    cat(sprintf("  individuals per population: %s\n",
                paste(range(x$n_ind, na.rm = TRUE), collapse = "-")))
  }
  cat("  populations:", paste(utils::head(x$pops, 8), collapse = ", "),
      if (length(x$pops) > 8) "..." else "", "\n")
  invisible(x)
}

# ---- internal helpers ---------------------------------------------------------

# Population assignment: a vector (one per individual, in file order), a
# data.frame or file with columns individual and population, or a function of
# the individual names.
resolve_pop <- function(pop, ind_names) {
  if (is.null(pop) || length(pop) == 0) {
    stop("Population assignments are needed. Supply `pop` as a vector with one ",
         "value per individual, a two-column data.frame or file (individual, ",
         "population), or a function that maps individual names to populations.",
         call. = FALSE)
  }
  if (is.function(pop)) pop <- pop(ind_names)
  if (is.character(pop) && length(pop) == 1 && length(ind_names) != 1 && file.exists(pop)) {
    pop <- utils::read.table(pop, header = FALSE, sep = "", colClasses = "character",
                             comment.char = "#", strip.white = TRUE)
  }
  if (is.data.frame(pop)) {
    map <- stats::setNames(as.character(pop[[2]]), as.character(pop[[1]]))
    if (!all(ind_names %in% names(map))) {
      if (all(map[1] == names(map)[1])) map <- map[-1]   # tolerate a header row
      miss <- setdiff(ind_names, names(map))
      if (length(miss)) stop(length(miss), " individuals have no population, e.g. ",
                             paste(utils::head(miss, 3), collapse = ", "), call. = FALSE)
    }
    pop <- map[ind_names]
  }
  if (length(pop) != length(ind_names)) {
    stop("`pop` has ", length(pop), " values but there are ", length(ind_names),
         " individuals.", call. = FALSE)
  }
  pop <- as.character(pop)
  factor(pop, levels = unique(pop))
}

# Allele calls (individual x locus x allele copy, character, NA for missing)
# to population allele counts
calls_to_popcounts <- function(calls, pop, loci) {
  counts <- lapply(seq_len(dim(calls)[2]), function(l) {
    a <- calls[, l, , drop = TRUE]
    if (is.null(dim(a))) a <- matrix(a, ncol = 1)
    al <- sort(unique(stats::na.omit(as.vector(a))))
    if (!length(al)) return(matrix(integer(0), nlevels(pop), 0))
    m <- matrix(0L, length(pop), length(al), dimnames = list(NULL, al))
    for (k in seq_len(ncol(a))) {
      ok <- !is.na(a[, k])
      idx <- cbind(which(ok), match(a[ok, k], al))
      m[idx] <- m[idx] + 1L
    }
    out <- rowsum(m, pop, reorder = FALSE)
    out[levels(pop), , drop = FALSE]
  })
  names(counts) <- loci
  new_popcounts(counts, levels(pop), table(pop)[levels(pop)])
}

# Allele dosage matrix (individual x SNP, 0..ploidy, NA missing) to counts
dosage_counts <- function(g, pop, ploidy = 2L) {
  alt <- rowsum(g, pop, reorder = FALSE, na.rm = TRUE)[levels(pop), , drop = FALSE]
  obs <- rowsum(ploidy * (!is.na(g)), pop, reorder = FALSE)[levels(pop), , drop = FALSE]
  loci <- if (is.null(colnames(g))) paste0("snp", seq_len(ncol(g))) else colnames(g)
  counts <- lapply(seq_len(ncol(g)), function(j) {
    cbind(ref = obs[, j] - alt[, j], alt = alt[, j])
  })
  names(counts) <- loci
  new_popcounts(counts, levels(pop), table(pop)[levels(pop)])
}

vcf_gt_counts <- function(calls, ids, pop) {
  calls[calls %in% c(".", "./.", ".|.")] <- NA
  split_gt <- strsplit(calls, "[/|]")
  ploidy <- max(lengths(split_gt), na.rm = TRUE)
  arr <- array(NA_character_, c(ncol(calls), nrow(calls), ploidy))
  for (k in seq_len(ploidy)) {
    v <- vapply(split_gt, function(z) if (length(z) >= k) z[k] else NA_character_, "")
    v[v == "."] <- NA
    arr[, , k] <- t(matrix(v, nrow(calls), ncol(calls)))
  }
  calls_to_popcounts(arr, pop, ids)
}
