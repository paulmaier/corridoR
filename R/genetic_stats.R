#' Pairwise FST between populations
#'
#' @param x A [popcounts][as_popcounts] object, or anything [as_popcounts()]
#'   accepts.
#' @param method `"hudson"` (default) is Hudson's FST as a ratio of averages
#'   across loci (Bhatia et al. 2013), which is robust to unequal sample sizes.
#'   `"nei"` is Nei's GST from pooled allele frequencies.
#' @return A symmetric matrix of pairwise FST with population names as
#'   dimnames. Any precomputed matrix (for example from \pkg{hierfstat} or
#'   Stacks) works just as well in the rest of the package.
#' @references Bhatia G, Patterson N, Sankararaman S, Price AL (2013)
#'   Estimating and interpreting FST: the impact of rare variants. Genome
#'   Research 23:1514-1521.
#' @examples
#' ex <- corridor_example()
#' fst <- pairwise_fst(read_genotypes(ex$files$small))
#' round(fst[1:4, 1:4], 3)
#' @export
pairwise_fst <- function(x, method = c("hudson", "nei")) {
  x <- as_popcounts(x)
  method <- match.arg(method)
  np <- length(x$pops)
  num <- den <- matrix(0, np, np)
  for (m in x$counts) {
    n <- rowSums(m)
    p <- m / n
    ok <- n > 1
    sq <- rowSums(p^2)
    if (method == "hudson") {
      hw <- (1 - sq) * n / (n - 1)                      # unbiased within-population diversity
      hb <- 1 - p %*% t(p)                              # between-population diversity
      w <- outer(hw, hw, "+") / 2
      a <- hb - w; b <- hb
    } else {
      hs <- outer(1 - sq, 1 - sq, "+") / 2
      pbar_sq <- (sq %o% rep(1, np) + rep(1, np) %o% sq + 2 * p %*% t(p)) / 4
      ht <- 1 - pbar_sq
      a <- ht - hs; b <- ht
    }
    use <- outer(ok, ok, "&")
    a[!use] <- 0; b[!use] <- 0
    num <- num + a; den <- den + b
  }
  out <- num / den
  diag(out) <- 0
  dimnames(out) <- list(x$pops, x$pops)
  out
}

#' Directional genetic differentiation and net migration (dM)
#'
#' Computes the directional GST of Sundqvist et al. (2016), as in
#' `diveRsity::divMigrate(stat = "gst")`: for each pair of populations, a
#' hypothetical shared gene pool is built from the normalized geometric mean of
#' their allele frequencies, and each population's differentiation from that
#' pool is converted to relative migration. Following Maier et al. (2022), the
#' difference between emigration and immigration gives net migration, dM.
#'
#' @inheritParams pairwise_fst
#' @return A list with
#'   \describe{
#'   \item{`relative`}{Relative migration, scaled so the largest value is 1.
#'     `relative[i, j]` is migration from population `i` into population `j`.}
#'   \item{`dM`}{`relative - t(relative)`: positive when `i` sends more migrants
#'     to `j` than it receives from it.}
#'   \item{`gst`}{The directional GST matrix behind `relative`.}
#'   }
#' @references Sundqvist L, Keenan K, Zackrisson M, Prodohl P, Kleinhans D
#'   (2016) Directional genetic differentiation and relative migration.
#'   Ecology and Evolution 6:3461-3475.
#'
#'   Maier PA, Vandergast AG, Ostoja SM, Aguilar A, Bohonak AJ (2022) Landscape
#'   genetics of a sub-alpine toad: climate change predicted to induce upward
#'   range shifts via asymmetrical migration corridors. Heredity 129:257-272.
#' @examples
#' ex <- corridor_example()
#' dg <- directional_gst(read_genotypes(ex$files$small))
#' round(dg$dM[1:4, 1:4], 3)
#' @export
directional_gst <- function(x) {
  x <- as_popcounts(x)
  np <- length(x$pops)
  nl <- length(x$counts)
  hs <- ht <- array(NA_real_, c(np, np, nl))
  for (l in seq_len(nl)) {
    m <- x$counts[[l]]
    n <- rowSums(m)
    af <- t(m / n)                                       # alleles x populations
    ok <- n > 0
    af[, !ok] <- 0
    # pooled gene pool f = sqrt(a * b) / sum(sqrt(a * b)), for all pairs at once
    q <- sqrt(af)
    S <- crossprod(q)                                    # sum_k sqrt(a_k b_k)
    f2 <- crossprod(af) / S^2                            # sum_k f_k^2
    fa <- crossprod(af^1.5, q) / S                       # sum_k f_k a_k, a from the row population
    a2 <- colSums(af^2)                                  # sum_k a_k^2, row population
    h_t <- 1 - (f2 + 2 * fa + a2) / 4
    h_s <- 1 - (f2 + a2) / 2
    h_t[!ok, ] <- NA; h_t[, !ok] <- NA; h_s[!ok, ] <- NA; h_s[, !ok] <- NA
    diag(h_t) <- NA; diag(h_s) <- NA
    ht[, , l] <- h_t; hs[, , l] <- h_s
  }
  mht <- apply(ht, c(1, 2), mean, na.rm = TRUE)
  mhs <- apply(hs, c(1, 2), mean, na.rm = TRUE)
  g <- (mht - mhs) / mht
  diag(g) <- 0
  mig <- (1 / g - 1) / 4
  mig[is.infinite(mig)] <- NA
  rel <- mig / max(mig, na.rm = TRUE)
  diag(rel) <- NA
  dn <- list(x$pops, x$pops)
  dimnames(g) <- dimnames(rel) <- dn
  dm <- rel - t(rel)
  diag(dm) <- 0
  list(relative = rel, dM = dm, gst = g)
}

#' Pairwise genetic table
#'
#' Turns matrices of FST and dM into a long table with one row per directed
#' pair, the layout used by the modeling functions.
#'
#' @param fst Symmetric matrix of pairwise differentiation.
#' @param dM Optional matrix of net migration, `dM[i, j]` from `i` to `j`.
#' @return A data.frame with columns `from`, `to`, `fst` and (if given) `dM`.
#' @examples
#' fst <- matrix(c(0, 0.1, 0.2, 0.1, 0, 0.15, 0.2, 0.15, 0), 3,
#'               dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
#' pairwise_table(fst)
#' @export
pairwise_table <- function(fst, dM = NULL) {
  ids <- rownames(fst)
  pr <- expand.grid(from = ids, to = ids, stringsAsFactors = FALSE)
  pr <- pr[pr$from != pr$to, ]
  out <- data.frame(from = pr$from, to = pr$to, fst = fst[cbind(pr$from, pr$to)])
  if (!is.null(dM)) out$dM <- dM[cbind(pr$from, pr$to)]
  rownames(out) <- NULL
  out
}
