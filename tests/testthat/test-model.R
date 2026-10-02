test_that("pair_weights down-weights frequent sites and sums to one", {
  w <- pair_weights(c("a", "a", "a", "b"), c("b", "c", "d", "c"))
  expect_equal(sum(w), 1)
  expect_lt(w[1], w[4])                      # pairs with the busy site "a" weigh less
})

test_that("vif_select keeps the more important of collinear features", {
  set.seed(1)
  x <- rnorm(100)
  X <- data.frame(a = x, b = x + rnorm(100, sd = 0.01), c = rnorm(100))
  expect_equal(sort(vif_select(X, c(a = 1, b = 2, c = 0.5), 10)), c("b", "c"))
})

test_that("select_bandwidth prefers the informative bandwidth", {
  skip_if_not_installed("ranger")
  set.seed(2)
  n <- 200
  good <- data.frame(from = "a", to = "b", snow = rnorm(n))
  bad <- data.frame(from = "a", to = "b", snow = rnorm(n))
  y <- good$snow + rnorm(n, sd = 0.2)
  sel <- select_bandwidth(list(bad = bad, good = good), y, list(climate = "snow"), num.trees = 100,
                          folds = 3, min.node.size = 5)
  expect_equal(unname(attr(sel, "best")), "good")
})

test_that("fit_connectivity fits, reports and projects", {
  skip_if_not_installed("Cubist")
  set.seed(3)
  n <- 250
  d <- data.frame(snow = rnorm(n), dist = runif(n))
  d$runoff <- d$snow + rnorm(n, sd = 0.2)
  d$temp <- -d$snow + rnorm(n, sd = 0.5)
  y <- 0.03 * d$snow + 0.05 * d$dist + rnorm(n, sd = 0.005)
  m <- fit_connectivity(d, y, list(climate = c("snow", "runoff", "temp")),
                        committees = c(1, 5), neighbors = c(0, 3), folds = 5)
  expect_s3_class(m, "corridor_model")
  expect_gt(m$best$Rsquared, 0.8)
  expect_true("dist" %in% m$vars)
  expect_output(print(m), "corridor_model")
  p_now <- predict(m, d)
  expect_length(p_now, n)
  p_fut <- predict(m, transform(d, snow = snow + 1, runoff = runoff + 1, temp = temp - 1))
  expect_gt(mean(p_fut - p_now), 0)          # more snow, higher response
})
