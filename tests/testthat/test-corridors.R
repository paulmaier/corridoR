# A small landscape with a wall across the middle and a gap at the top
small_world <- function() {
  crs <- "+proj=tmerc +lat_0=0 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  r <- terra::rast(nrows = 30, ncols = 40, xmin = 0, xmax = 4000, ymin = 0, ymax = 3000, crs = crs)
  terra::values(r) <- 1
  wall <- terra::colFromX(r, 2000)
  r[5:30, wall] <- 1e6
  sites <- sf::st_as_sf(data.frame(site = c("w", "e", "s"), x = c(500, 3500, 2500), y = c(2500, 2500, 500)),
                        coords = c("x", "y"), crs = crs)
  list(r = r, sites = sites)
}

test_that("paths avoid barriers and corridors are weighted 0 to 1", {
  w <- small_world()
  tr <- make_transition(w$r, barrier = 1e6)
  p <- least_cost_paths(tr, w$sites, pairs = data.frame(from = "w", to = "s"))
  expect_gt(p$length, p$euclidean)                     # forced through the gap
  acc <- accumulated_cost(tr, w$sites)
  expect_equal(names(acc), c("w", "e", "s"))
  lw <- lcc_weights(acc, "w", "e", q = c(0.01, 0.1))
  v1 <- terra::values(lw[[1]], na.rm = TRUE); v2 <- terra::values(lw[[2]], na.rm = TRUE)
  expect_true(all(v1 >= 0 & v1 <= 1) && all(v2 >= 0 & v2 <= 1))
  expect_lt(length(v1), length(v2))                    # the broader corridor keeps more cells
  expect_equal(max(v2), 1)
})

test_that("corridor_extract gives weighted means and sums", {
  w <- small_world()
  tr <- make_transition(w$r, barrier = 1e6)
  acc <- accumulated_cost(tr, w$sites)
  env <- c(w$r * 0 + 7, w$r * 0 + 1)
  names(env) <- c("constant", "count")
  pr <- data.frame(from = c("w", "w"), to = c("e", "s"))
  x <- corridor_extract(env, pr, acc = acc, q = 0.05, stat = c(count = "sum"))
  expect_equal(x$constant, c(7, 7))
  wts <- terra::values(lcc_weights(acc, "w", "e", 0.05), na.rm = TRUE)
  expect_equal(x$count[1], sum(wts))
  p <- least_cost_paths(tr, w$sites, pairs = pr)
  b <- path_buffers(p, 200)[[1]]
  xb <- corridor_extract(env, pr, buffers = b)
  expect_equal(xb$constant, c(7, 7))
})

test_that("site_contrast is source minus destination", {
  w <- small_world()
  env <- terra::init(w$r, "x")
  names(env) <- "east"
  pr <- both_directions(data.frame(from = "w", to = "e"))
  sc <- site_contrast(w$sites, pr, env = env, attributes = data.frame(site = c("w", "e", "s"), area = c(1, 5, 2)))
  expect_equal(sc$east.at, c(-3000, 3000))
  expect_equal(sc$area.at, c(-4, 4))
})

test_that("map_pairwise points arrows from source to destination", {
  w <- small_world()
  tr <- make_transition(w$r, barrier = 1e6)
  acc <- accumulated_cost(tr, w$sites)
  v <- data.frame(from = c("w", "e"), to = c("e", "w"), value = c(0.2, -0.2))
  m <- map_pairwise(v, acc, w$sites)
  expect_equal(names(m), c("value", "vx", "vy", "bearing"))
  b <- terra::values(m$bearing, na.rm = TRUE)
  expect_true(all(abs(b - 90) < 1e-6))                 # due east
  expect_equal(max(terra::values(m$value, na.rm = TRUE)), 0.2)
  a <- shift_arrows(m, cells = 5, min_quantile = 0)
  expect_true(all(a$xend > a$x))
})

test_that("plot_shift draws without error", {
  w <- small_world()
  tr <- make_transition(w$r, barrier = 1e6)
  acc <- accumulated_cost(tr, w$sites)
  m <- map_pairwise(data.frame(from = c("w", "e"), to = c("e", "w"), value = c(0.2, -0.2)), acc, w$sites)
  f <- tempfile(fileext = ".png"); grDevices::png(f)
  expect_invisible(plot_shift(m, dem = w$r, sites = w$sites, cells = 5))
  grDevices::dev.off()
})
