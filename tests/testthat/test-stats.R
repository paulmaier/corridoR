counts <- function(...) {
  m <- rbind(...)
  rownames(m) <- paste0("P", seq_len(nrow(m)))
  m
}
mk <- function(loci) new_popcounts(loci, rownames(loci[[1]]), rep(10, nrow(loci[[1]])))

test_that("Hudson's FST matches a hand calculation", {
  pc <- mk(list(counts(c(15, 5), c(5, 15))))
  n <- 20; p1 <- 0.75; p2 <- 0.25
  hw <- 2 * p1 * (1 - p1) * n / (n - 1)
  hb <- p1 * (1 - p2) + p2 * (1 - p1)
  expect_equal(pairwise_fst(pc)[1, 2], 1 - hw / hb)
  expect_equal(pairwise_fst(pc), t(pairwise_fst(pc)))
})

test_that("identical populations have FST near zero", {
  pc <- mk(list(counts(c(10, 10), c(10, 10)), counts(c(4, 16), c(4, 16))))
  expect_lt(pairwise_fst(pc)[1, 2], 0)                # unbiased estimator dips below zero
  expect_equal(pairwise_fst(pc, "nei")[1, 2], 0)
})

test_that("directional GST follows the divMigrate formula", {
  pc <- mk(list(counts(c(12, 8), c(4, 16), c(10, 10)), counts(c(20, 0), c(14, 6), c(9, 11))))
  dg <- directional_gst(pc)
  # hand calculation for P1 relative to the pool with P2, averaged over loci
  part <- function(a, b) {
    f <- sqrt(a * b) / sum(sqrt(a * b))
    c(ht = 1 - sum(((f + a) / 2)^2), hs = 1 - sum((f^2 + a^2) / 2))
  }
  af <- lapply(pc$counts, function(m) m / rowSums(m))
  h <- rowMeans(sapply(af, function(x) part(x[1, ], x[2, ])))
  g12 <- (h["ht"] - h["hs"]) / h["ht"]
  expect_equal(unname(dg$gst[1, 2]), unname(g12))
  mig <- (1 / dg$gst - 1) / 4
  diag(mig) <- NA
  expect_equal(dg$relative, mig / max(mig, na.rm = TRUE))
  expect_equal(dg$dM, -t(dg$dM))
  expect_equal(max(dg$relative, na.rm = TRUE), 1)
})

test_that("directional GST recovers the direction of simulated migration", {
  ex <- corridor_example()
  pc <- read_genotypes(ex$files$structure)
  wet <- terra::extract(ex$env$moisture, terra::vect(ex$sites))[, 2]
  names(wet) <- ex$sites$site
  tab <- pairwise_table(pairwise_fst(pc), directional_gst(pc)$dM)
  # migrants favor wetter meadows, so dM rises with the destination's moisture
  expect_gt(stats::cor(tab$dM, wet[tab$to] - wet[tab$from]), 0.3)
})

test_that("pairwise_table lists every directed pair", {
  m <- matrix(c(0, .1, .2, .1, 0, .3, .2, .3, 0), 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  t <- pairwise_table(m, m - t(m))
  expect_equal(nrow(t), 6)
  expect_equal(t$fst[t$from == "a" & t$to == "c"], 0.2)
})
