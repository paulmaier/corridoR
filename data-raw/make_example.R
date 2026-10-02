# Builds the simulated example shipped in inst/extdata.
#
# A 60 x 40 km range modeled loosely on the central Sierra Nevada: a long,
# gradual western slope rising from foothills to a crest near the eastern
# edge, a steep eastern escarpment in a rain shadow, and two west-flowing
# river canyons that deepen downstream. Eighty meadows sit between about 1000
# and 3400 m, a few of them on benches below the crest on the eastern side.
# Four lineages, arranged as in Yosemite: two at low elevation on either side
# of the southern canyon (West and South) and two at high elevation (North,
# above the northern canyon, and East along the crest). Lineages at the same
# elevation share their environment but not their history.
#
# Toads move most easily over gentle, forested ground. Migrants favor meadows
# that are wetter *relative to* the one they leave (the moisture index is on a
# log scale), so wetter, higher meadows receive more migrants than they send.
# Under warming (+4.2 C) the snow line climbs about 500 m. Low meadows lose the
# largest share of their moisture and the crest the smallest, so the pull
# toward higher meadows strengthens most in the foothills and fades toward the
# crest. Each lineage has its own moisture requirement: how strongly its
# migrants respond to a wetter destination differs between lineages.
#
# Run from the package root: Rscript data-raw/make_example.R

suppressPackageStartupMessages({ library(terra); library(sf) })
devtools::load_all(quiet = TRUE)
set.seed(2027)
out <- "inst/extdata"
crs_local <- "+proj=tmerc +lat_0=0 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

# ---- terrain ------------------------------------------------------------------
r <- rast(nrows = 200, ncols = 300, xmin = 0, xmax = 60000, ymin = 0, ymax = 40000, crs = crs_local)
xy <- xyFromCell(r, seq_len(ncell(r)))
x <- xy[, 1] / 60000; y <- xy[, 2] / 40000

smooth_noise <- function(w, sd = 1) {
  n <- setValues(r, rnorm(ncell(r)))
  n <- focal(n, w = focalMat(n, w, "Gauss"), fun = "sum", na.policy = "omit", fillvalue = 0)
  n / global(n, "sd", na.rm = TRUE)[[1]] * sd
}
crest <- 0.82 + 0.025 * sin(2 * pi * 1.4 * y + 1)            # the crest wanders a little
west <- x <= crest
base <- ifelse(west,
               1150 + 2450 * (x / crest)^1.7,                 # long gradual western slope
               3600 - 2000 * ((x - crest) / (1 - crest))^0.45) # steep eastern escarpment
spurs <- 220 * sin(2 * pi * (3.5 * y + 0.8 * x)) * pmin(1, x / crest)
canyon <- function(y0, depth, head = 0.72) {
  down <- pmin(1, pmax(0, (head - x) / head))^0.6              # deepens downstream (west)
  wid <- 0.008 + 0.012 * down                                   # narrow, steep-walled gorge
  yc <- y0 + 0.035 * sin(8 * x + y0 * 10)
  -depth * down * exp(-((y - yc)^2) / (2 * wid^2))
}
elev <- base + spurs + canyon(0.33, 1300) + canyon(0.70, 1150)
dem <- setValues(r, elev) + smooth_noise(1800, 110) + smooth_noise(600, 35)
names(dem) <- "elevation"
slope <- terrain(dem, "slope", unit = "degrees")
slope <- focal(slope, 3, mean, na.policy = "only", na.rm = TRUE)

# ---- climate, present and future ------------------------------------------------
logistic <- function(z) 1 / (1 + exp(-z))
shadow <- setValues(r, ifelse(west, 1, 0.55))                 # rain shadow east of the crest
shadow <- focal(shadow, focalMat(shadow, 1200, "Gauss"), fun = "sum", na.policy = "omit", fillvalue = 1)
# April 1 snowpack rises with elevation. Warming removes a share of it that
# falls smoothly with elevation, as in the Basin Characterization Model for
# Yosemite (about 90% near 1400 m, under 20% above 3100 m).
n_snow <- smooth_noise(2500, 0.06); n_run <- smooth_noise(3000, 10)
snow_now <- 1600 * logistic((dem - 1700) / 520) * shadow * exp(n_snow)
loss <- clamp(0.95 * exp(-(dem - 1350) / 950), 0, 0.95)
snow_fut <- snow_now * (1 - loss)
runoff_now <- clamp(0.6 * snow_now + 60 * shadow + n_run, 1)
runoff_fut <- clamp(0.6 * snow_fut + 60 * shadow + n_run, 1)
# meadow moisture index on a log scale: a change is relative to how wet a
# meadow is now, so the meadows losing the largest share change the most
moisture_now <- log(snow_now + runoff_now + 20)
moisture_fut <- log(snow_fut + runoff_fut + 20)
temp_now <- 33 - 6.4 * dem / 1000 + 1.5 * (1 - shadow) + smooth_noise(4000, 0.4)
temp_fut <- temp_now + 4.2
forest <- clamp(0.9 * logistic(-(dem - 2800) / 260) * (0.6 + 0.4 * shadow) + smooth_noise(1200, 0.12), 0, 1)

env <- c(snow_now, runoff_now, moisture_now, temp_now, slope, forest)
names(env) <- c("snowpack", "runoff", "moisture", "summer_temp", "slope", "forest")
env_future <- c(snow_fut, runoff_fut, moisture_fut, temp_fut, slope, forest)
names(env_future) <- names(env)

# ---- resistance: two components, blended 80:20 in the "true" model -------------------
res_slope <- 1 + 99 * (1 - exp(-slope / 16))^2
res_cover <- 1 + 99 * (1 - forest)^1.5
resistance <- 0.8 * res_slope + 0.2 * res_cover
# the escarpment is impassable except at three passes over the crest
passes <- c(0.18, 0.5, 0.83)
near_pass <- setValues(r, apply(abs(outer(y, passes, "-")), 1, min) < 0.03)
ridge <- slope > 33 & !near_pass
resistance[ridge] <- 1e6; res_slope[ridge] <- 1e6; res_cover[ridge] <- 1e6
names(resistance) <- "resistance"; names(res_slope) <- "resistance_slope"; names(res_cover) <- "resistance_cover"

# ---- meadows: 14 in the foothills, 52 on the western slope, 14 east of the crest ----
ok <- values(slope)[, 1] < 9 & values(dem)[, 1] < 3450 & y > 0.05 & y < 0.95 & x > 0.04 & x < 0.97
el <- values(dem)[, 1]
place <- function(cand, n, gap, have = integer(0)) {
  pick <- have
  for (c0 in sample(cand)) {
    if (!length(pick) || min(sqrt(colSums((t(xy[pick, , drop = FALSE]) - xy[c0, ])^2))) > gap) pick <- c(pick, c0)
    if (length(pick) == length(have) + n) break
  }
  pick
}
pick <- place(which(ok & !west & el > 1700), 14, 2200)
pick <- place(which(ok & west & el > 1000 & el < 1500 & x < 0.3), 14, 3000, pick)
pick <- place(which(ok & west & el >= 1500), 52, 3000, pick)
pick <- pick[order(xy[pick, 1])]
sites <- st_as_sf(data.frame(site = sprintf("M%02d", seq_along(pick)), xy[pick, , drop = FALSE]),
                  coords = c("x", "y"), crs = crs_local)
sites$elevation <- round(extract(dem, vect(sites))[, 2])
sxy <- st_coordinates(sites); sx <- sxy[, 1] / 60000; sy <- sxy[, 2] / 40000
canyon_y <- function(y0, x) y0 + 0.035 * sin(8 * x + y0 * 10)          # canyon centerlines
high <- sites$elevation >= 2350 | sx > 0.82 + 0.025 * sin(2 * pi * 1.4 * sy + 1)
sites$lineage <- ifelse(high, ifelse(sy > canyon_y(0.70, sx), "North", "East"),
                        ifelse(sy > canyon_y(0.33, sx), "West", "South"))
# A meadow surrounded by another lineage looks odd on a map; swap the most
# out-of-place South and East meadows so each lineage forms one cluster
nn_share <- function(i, l) {
  d <- sqrt(colSums((t(sxy) - sxy[i, ])^2)); nb <- order(d)[2:6]
  mean(sites$lineage[nb] == l)
}
odd_s <- which(sites$lineage == "South"); odd_e <- which(sites$lineage == "East")
s_out <- odd_s[which.max(sapply(odd_s, nn_share, l = "East"))]
e_out <- odd_e[which.max(sapply(odd_e, nn_share, l = "South"))]
if (nn_share(s_out, "East") > 0.5) {
  sites$lineage[c(s_out, e_out)] <- sites$lineage[c(e_out, s_out)]
}
message(nrow(sites), " meadows, ", paste(range(sites$elevation), collapse = "-"), " m; lineages: ",
        paste(names(table(sites$lineage)), table(sites$lineage), collapse = ", "))

# ---- migration matrix ----------------------------------------------------------------
tr <- make_transition(resistance, barrier = 1e6)
acc <- accumulated_cost(tr, sites)
cost <- sapply(sites$site, function(i) extract(acc[[i]], vect(sites))[, 2])
dimnames(cost) <- list(sites$site, sites$site)
wet <- extract(moisture_now, vect(sites))[, 2]
wet_z <- (wet - mean(wet)) / sd(wet)
np <- nrow(sites)
lin <- sites$lineage
# Water stress: how far a meadow falls short of the wettest meadow today.
# Migrants leave a stressed meadow for wetter neighbors; the more stressed
# their meadow, the stronger the pull. Each lineage has its own strength.
beta <- c(West = 1.3, South = 1.4, North = 1.1, East = 1.2)[lin]
wmax <- max(wet)
stress <- function(w) pmax(0, (wmax - w) / sd(wet))
base <- exp(-cost / quantile(cost[is.finite(cost) & cost > 0], 0.12))  # m[i, j]: into i from j
base[!is.finite(base)] <- 0
base[outer(lin, lin, "!=")] <- base[outer(lin, lin, "!=")] * 0.9      # slightly fewer between lineages
pull <- function(w) {
  dz <- outer(w, w, "-") / sd(wet)                                     # destination i wetter than source j
  # the pull saturates, so very dry meadows cannot send unlimited migrants
  base * exp(2 * tanh(sweep(pmax(dz, 0), 2, beta * stress(w), "*") / 3) - 0.3 * pmax(-dz, 0))
}
raw_now <- pull(wet)
diag(raw_now) <- 0
k_mig <- 0.06 / max(rowSums(raw_now))
m <- raw_now * k_mig
diag(m) <- 1 - rowSums(m)

# ---- Wright-Fisher drift and migration, 1000 biallelic SNPs -----------------------------
nloc <- 1500; N <- 150; gens <- 500
# lineages start from diverged ancestral gene pools (Balding-Nichols, F = 0.02)
p0 <- rbeta(nloc, 0.6, 0.6)
anc <- sapply(unique(lin), function(l) rbeta(nloc, p0 * (1 / 0.02 - 1), (1 - p0) * (1 / 0.02 - 1)))
p <- t(anc[, match(lin, unique(lin))])
for (g in seq_len(gens)) {
  p <- m %*% p
  p[] <- rbinom(length(p), 2 * N, p) / (2 * N)
}
n_ind <- 8
geno <- do.call(rbind, lapply(seq_len(np), function(i) {
  matrix(rbinom(n_ind * nloc, 2, rep(p[i, ], each = n_ind)), n_ind, nloc)
}))
geno <- geno[, apply(geno, 2, function(v) length(unique(v)) > 1)]
geno[sample(length(geno), round(0.03 * length(geno)))] <- NA          # 3% missing calls
ind <- sprintf("%s_%02d", rep(sites$site, each = n_ind), rep(seq_len(n_ind), np))
pop <- rep(sites$site, each = n_ind)
loci <- sprintf("snp%03d", seq_len(ncol(geno)))
message(ncol(geno), " polymorphic SNPs for ", length(ind), " individuals")

# ---- genotype files ---------------------------------------------------------------------
dir.create(out, showWarnings = FALSE, recursive = TRUE)
a1 <- ifelse(is.na(geno), -9, ifelse(geno == 2, 2, 1))     # STRUCTURE: alleles coded 1 and 2
a2 <- ifelse(is.na(geno), -9, ifelse(geno >= 1, 2, 1))
writeLines(c(paste(loci, collapse = " "),
             unlist(lapply(seq_along(ind), function(k) c(
               paste(ind[k], pop[k], paste(a1[k, ], collapse = " ")),
               paste(ind[k], pop[k], paste(a2[k, ], collapse = " ")))))),
           file.path(out, "example.str"))
gp <- ifelse(is.na(geno), "0000", ifelse(geno == 0, "0101", ifelse(geno == 1, "0102", "0202")))
gen_lines <- c("Simulated toad meadows (corridoR example)", loci)
for (s in sites$site) {
  k <- which(pop == s)
  gen_lines <- c(gen_lines, "Pop", paste0(ind[k], " , ", apply(gp[k, , drop = FALSE], 1, paste, collapse = " ")))
}
writeLines(gen_lines, file.path(out, "example.gen"))
gt <- ifelse(is.na(geno), "./.", ifelse(geno == 0, "0/0", ifelse(geno == 1, "0/1", "1/1")))
writeLines(c("##fileformat=VCFv4.2", "##source=corridoR simulated example",
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", ind), collapse = "\t"),
             vapply(seq_along(loci), function(j) paste(c("chr1", j * 1000, loci[j], "A", "G", ".", "PASS", ".", "GT",
                                                         gt[, j]), collapse = "\t"), "")),
           file.path(out, "example.vcf"))
write.table(data.frame(ind, pop), file.path(out, "popmap.txt"), quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = "\t")

# ---- rasters and sites -------------------------------------------------------------------
gd <- c("COMPRESS=DEFLATE", "PREDICTOR=3")
writeRaster(c(dem, resistance, res_slope, res_cover), file.path(out, "landscape.tif"), overwrite = TRUE,
            gdal = gd, datatype = "FLT4S")
writeRaster(env, file.path(out, "env_present.tif"), overwrite = TRUE, gdal = gd, datatype = "FLT4S")
writeRaster(env_future, file.path(out, "env_future.tif"), overwrite = TRUE, gdal = gd, datatype = "FLT4S")
st_write(sites, file.path(out, "sites.gpkg"), delete_dsn = TRUE, quiet = TRUE)
# time since each pair of lineages split (thousand years), as in Maier et al. (2019)
write.csv(data.frame(lineage1 = c("North", "North", "North", "East", "East", "West"),
                     lineage2 = c("East", "West", "South", "West", "South", "South"),
                     tmrca = c(720, 720, 720, 300, 420, 420)),
          file.path(out, "lineage_tmrca.csv"), row.names = FALSE)
# the same migration model under future moisture, for checking forecasts
wet_f <- extract(moisture_fut, vect(sites))[, 2]
raw_fut <- pull(wet_f); diag(raw_fut) <- 0
saveRDS(list(migration_now = raw_now * k_mig, migration_future = raw_fut * k_mig, wet = wet,
             wet_future = wet_f, beta = beta), "data-raw/simulation_truth.rds")
message("done")
