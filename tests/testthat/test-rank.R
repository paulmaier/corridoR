test_that("rank_resistance puts the surface that shaped genetic distance first", {
  skip_if_not_installed("lme4")
  set.seed(4)
  r <- terra::rast(nrows = 30, ncols = 30, xmin = 0, xmax = 3000, ymin = 0, ymax = 3000,
                   crs = "+proj=utm +zone=11 +datum=WGS84 +units=m +no_defs")
  good <- terra::setValues(r, rep(c(1, 1, 1, 50, 50), length.out = 900))   # bands of high cost
  bad <- terra::setValues(r, runif(900, 1, 50))
  xy <- cbind(runif(12, 100, 2900), runif(12, 100, 2900))
  sites <- sf::st_as_sf(data.frame(site = sprintf("s%02d", 1:12), x = xy[, 1], y = xy[, 2]),
                        coords = c("x", "y"), crs = terra::crs(r))
  cd <- as.matrix(gdistance::costDistance(make_transition(good), xy))
  gen <- 0.02 + cd / max(cd) * 0.1 + matrix(rnorm(144, sd = 0.002), 12)
  gen <- (gen + t(gen)) / 2
  dimnames(gen) <- list(sites$site, sites$site)
  res <- rank_resistance(list(good = good, bad = bad), sites, gen)
  expect_equal(res$hypothesis[1], "good")
  expect_setequal(res$hypothesis, c("good", "bad", "straight_line"))
  expect_equal(res$delta_AIC[1], 0)
  expect_true(all(res$R2m >= 0 & res$R2m <= 1))
  expect_equal(nrow(attr(res, "distances")), choose(12, 2))
  # a multi-layer SpatRaster works the same way
  both <- c(good, bad); names(both) <- c("good", "bad")
  expect_equal(rank_resistance(both, sites, gen, straight = FALSE)$hypothesis[1], "good")
})
