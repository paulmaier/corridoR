# One validation replicate: simulate a landscape, run the corridoR workflow,
# and compare forecasts of change in FST and dM with the change built into the
# simulation. Four alternatives are fit for comparison: Cubist with every
# feature group in a 400 m path buffer, Cubist with straight-line transects,
# Cubist with every group in the same least cost corridor (q = 0.01), and a
# random forest on the corridoR features.
#
# Usage (from the package root):
#   Rscript data-raw/validation/run_replicate.R <seed> <outdir> <scenario>
# with scenario "valleys" or "valleys_warmer" (see data-raw/make_example.R).
#
# Settings (environment variables):
#   REUSE=1        analyse a landscape already simulated in <outdir>
#   RESIST         surface the analysis uses: "true" (the one that generated
#                  gene flow, default) or "slope" (its slope component only,
#                  which misses that toads also track moisture)
#   PATCHES        patchy variation added to the moisture layer on the land
#                  between meadows, in units of the layer's standard deviation
#                  (default 0)
#   PATCH_SIZE     size of those patches in metres (default 1500)
args <- commandArgs(TRUE)
seed <- as.integer(args[1]); dir <- args[2]; scen <- args[3]
dir.create(dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(SIM_SEED = seed, SIM_OUT = dir, SIM_SCENARIO = scen)
resist <- Sys.getenv("RESIST", "true")
patches <- as.numeric(Sys.getenv("PATCHES", "0"))
patch_size <- as.numeric(Sys.getenv("PATCH_SIZE", "1500"))
tag <- paste0("_", resist, if (patches > 0) sprintf("_patches%g_%gm", patches, patch_size) else "")
if (Sys.getenv("REUSE") == "1" && file.exists(file.path(dir, "simulation_truth.rds"))) {
  devtools::load_all(quiet = TRUE)
} else {
  source("data-raw/make_example.R")
}
suppressPackageStartupMessages({ library(terra); library(sf) })
t0 <- proc.time()

f <- function(x) file.path(dir, x)
land <- rast(f("landscape.tif"))
dem <- land[["elevation"]]
resistance <- land[[c(true = "resistance", slope = "resistance_slope")[[resist]]]]
env <- rast(f("env_present.tif")); env_future <- rast(f("env_future.tif"))
# Patches stand for rock, water and canopy in a satellite moisture index: the
# same today and in future, and absent from the true moisture that shapes gene
# flow. Meadows themselves are measured from their own, consistently wet
# polygons, so the patches affect only the route measurements.
env_site <- env; env_future_site <- env_future
if (patches > 0) {
  set.seed(seed + 9999)
  nz <- env$moisture; values(nz) <- rnorm(ncell(nz))
  nz <- focal(nz, focalMat(nz, patch_size, "Gauss"), fun = "sum", na.policy = "omit", fillvalue = 0)
  nz <- nz / global(nz, "sd", na.rm = TRUE)[[1]] * patches * global(env$moisture, "sd", na.rm = TRUE)[[1]]
  env$moisture <- env$moisture + nz; env_future$moisture <- env_future$moisture + nz
}
sites <- st_read(f("sites.gpkg"), quiet = TRUE)
tmrca <- read.csv(f("lineage_tmrca.csv"))
truth <- readRDS(f("simulation_truth.rds"))

# Step 1: genetics
pc <- read_genotypes(f("example.str"))
fst <- pairwise_fst(pc); gen <- pairwise_table(fst, directional_gst(pc)$dM)

# Steps 2-3: paths, corridors, extraction
tr <- make_transition(resistance, barrier = 1e6)
acc <- accumulated_cost(tr, sites)
xy <- st_coordinates(sites); rownames(xy) <- sites$site
pairs <- as.data.frame(t(combn(sites$site, 2))); names(pairs) <- c("from", "to")
pairs <- pairs[sqrt(rowSums((xy[pairs$from, ] - xy[pairs$to, ])^2)) < 15000, ]
paths <- least_cost_paths(tr, sites, pairs = pairs, dem = dem)
both <- c(env, env_future); names(both) <- c(names(env), paste0(names(env), "_f"))
# straight-line transects between meadows, the usual alternative to paths
seg <- st_sfc(lapply(seq_len(nrow(pairs)), function(k)
  st_linestring(rbind(xy[pairs$from[k], ], xy[pairs$to[k], ]))), crs = st_crs(sites))
transects <- st_sf(from = pairs$from, to = pairs$to, geometry = st_buffer(seg, 400))
transect <- corridor_extract(both, pairs, buffers = transects, progress = FALSE)
bands <- list(
  lcp_400m = corridor_extract(both, pairs, buffers = path_buffers(paths, 400)[[1]], progress = FALSE),
  lcc_0.01 = corridor_extract(both, pairs, acc = acc, q = 0.01, progress = FALSE),
  lcc_0.05 = corridor_extract(both, pairs, acc = acc, q = 0.05, progress = FALSE))

# Step 4: bandwidth per group
groups <- list(climate = c("snowpack", "runoff", "summer_temp"), wetness = "moisture",
               terrain = "slope", cover = "forest")
sel <- select_bandwidth(bands, fst[cbind(pairs$from, pairs$to)], groups, num.trees = 300)
best <- attr(sel, "best")
write.csv(sel, f(paste0("bandwidth", tag, ".csv")), row.names = FALSE)

build <- function(bw_of, straight = FALSE) {
  pick <- function(suffix = "") {
    out <- pairs
    src <- if (straight) list(t = transect) else bands
    for (g in names(groups)) out[groups[[g]]] <- src[[if (straight) "t" else bw_of[[g]]]][paste0(groups[[g]], suffix)]
    out$path_length <- if (straight) sqrt(rowSums((xy[pairs$from, ] - xy[pairs$to, ])^2)) else paths$length
    both_directions(out)
  }
  now <- pick(); fut <- pick("_f")
  now <- cbind(now, site_contrast(sites, now, env = env_site)[, -(1:2)])
  fut <- cbind(fut, site_contrast(sites, now, env = env_future_site)[, -(1:2)])
  lin <- st_drop_geometry(sites)[, c("site", "lineage")]
  now$LineageCross <- fut$LineageCross <- lineage_cross(now$from, now$to, lin, tmrca)
  list(now = now, fut = fut)
}
d <- build(best)
d_buf <- build(setNames(rep("lcp_400m", length(groups)), names(groups)))
d_tr <- build(NULL, straight = TRUE)
d_lcc <- build(setNames(rep("lcc_0.01", length(groups)), names(groups)))   # corridors for every group
y <- gen$dM[match(paste(d$now$from, d$now$to), paste(gen$from, gen$to))]
mg <- c(groups, list(climate.at = paste0(groups$climate, ".at")))

# Step 5-6: models and forecasts, for dM (direction) and FST (amount of gene flow)
yf <- gen$fst[match(paste(d$now$from, d$now$to), paste(gen$from, gen$to))]
feat <- setdiff(names(d$now), c("from", "to"))
forecast <- function(resp) {
  m <- fit_connectivity(d$now, resp, mg, committees = c(1, 10, 50), neighbors = c(0, 5))
  m_buf <- fit_connectivity(d_buf$now, resp, mg, committees = c(1, 10, 50), neighbors = c(0, 5))
  m_tr <- fit_connectivity(d_tr$now, resp, mg, committees = c(1, 10, 50), neighbors = c(0, 5))
  m_lcc <- fit_connectivity(d_lcc$now, resp, mg, committees = c(1, 10, 50), neighbors = c(0, 5))
  rf <- ranger::ranger(x = d$now[, feat], y = resp, num.trees = 1000, seed = seed, num.threads = 2)
  list(change = data.frame(
         corridoR = predict(m, d$fut) - predict(m, d$now),
         buffer_only = predict(m_buf, d_buf$fut) - predict(m_buf, d_buf$now),
         transect = predict(m_tr, d_tr$fut) - predict(m_tr, d_tr$now),
         corridor_fixed = predict(m_lcc, d_lcc$fut) - predict(m_lcc, d_lcc$now),
         random_forest = predict(rf, d$fut[, feat])$predictions - predict(rf, d$now[, feat])$predictions),
       cv = c(m$best$Rsquared, m_buf$best$Rsquared, m_tr$best$Rsquared, m_lcc$best$Rsquared, NA))
}
fd <- forecast(y); ff <- forecast(yf)

# Truth: change in net migration (for dM) and in total gene flow (for FST;
# more gene flow means lower FST, so its sign is flipped)
a <- d$now$from; b <- d$now$to
m0 <- truth$migration_now; m1 <- truth$migration_future
# dM describes asymmetry, so its truth is the change in the share of a pair's
# flow going each way; FST describes the amount of flow (as in Maier et al. 2022,
# which models magnitude and asymmetry separately). The change in absolute net
# flow, which mixes the two, is also reported.
share <- function(m) (m[cbind(b, a)] - m[cbind(a, b)]) / (m[cbind(b, a)] + m[cbind(a, b)])
true_dm <- share(m1) - share(m0)
true_dm_abs <- (m1[cbind(b, a)] - m1[cbind(a, b)]) - (m0[cbind(b, a)] - m0[cbind(a, b)])
true_fst <- -((m1[cbind(b, a)] + m1[cbind(a, b)]) - (m0[cbind(b, a)] + m0[cbind(a, b)]))
score <- function(pred, truth) {
  ok <- is.finite(pred) & is.finite(truth); pred <- pred[ok]; truth <- truth[ok]
  big <- abs(truth) >= quantile(abs(truth), 0.75)
  c(r = cor(pred, truth), sign_top25 = mean(sign(pred[big]) == sign(truth[big])))
}
methods <- c("corridoR", "buffer_only", "transect", "corridor_fixed", "random_forest")
res <- rbind(
  data.frame(response = "dM", method = methods,
             t(sapply(methods, function(k) score(fd$change[[k]], true_dm))), cv_r2 = fd$cv),
  data.frame(response = "dM_absolute", method = methods,
             t(sapply(methods, function(k) score(fd$change[[k]], true_dm_abs))), cv_r2 = fd$cv),
  data.frame(response = "FST", method = methods,
             t(sapply(methods, function(k) score(ff$change[[k]], true_fst))), cv_r2 = ff$cv))
res <- cbind(scenario = scen, seed = seed, outside_range = truth$outside_range,
             best_climate_bw = best[["climate"]], best_wetness_bw = best[["wetness"]],
             best_terrain_bw = best[["terrain"]], best_cover_bw = best[["cover"]], res,
             dM_vs_truth_net = cor(gen$dM, gen_net <- (m0[cbind(gen$to, gen$from)] - m0[cbind(gen$from, gen$to)])),
             minutes = (proc.time() - t0)[3] / 60, row.names = NULL)
res$resist <- resist; res$patches <- patches; res$patch_size <- patch_size
write.csv(res, f(paste0("result", tag, ".csv")), row.names = FALSE)
saveRDS(list(from = a, to = b, dM = fd$change, FST = ff$change, true_dm = true_dm, true_dm_abs = true_dm_abs, true_fst = true_fst), f(paste0("change", tag, ".rds")))
print(res)
