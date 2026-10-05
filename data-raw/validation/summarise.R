# Summarises the validation study (see run_validation.sh): medians and
# 10th-90th percentiles of forecast accuracy, and paired Wilcoxon tests of
# corridoR against each alternative across landscapes.
# Usage: Rscript data-raw/validation/summarise.R <output folder>
root <- commandArgs(TRUE)[1]
read_set <- function(scen, file, label) {
  f <- list.files(file.path(root, scen), paste0("^", file, "$"), recursive = TRUE, full.names = TRUE)
  if (!length(f)) return(NULL)
  x <- do.call(rbind, lapply(f, read.csv))
  x$set <- label
  x
}
x <- rbind(read_set("valleys", "result_true.csv", "A"),
           read_set("valleys", "result_slope_patches1_1500m.csv", "B"),
           read_set("valleys_warmer", "result_slope_patches1_1500m.csv", "C"))
x <- x[x$response %in% c("FST", "dM"), ]
methods <- c("corridoR", "buffer_only", "transect", "corridor_fixed", "random_forest")
bw <- unique(x[, c("set", "seed", "best_climate_bw", "best_wetness_bw", "best_terrain_bw", "best_cover_bw")])
long <- do.call(rbind, lapply(split(x, list(x$set, x$response, x$method), drop = TRUE), function(d) {
  ref <- x$r[x$set == d$set[1] & x$response == d$response[1] & x$method == "corridoR"][order(x$seed[x$set == d$set[1] & x$response == d$response[1] & x$method == "corridoR"])]
  v <- d$r[order(d$seed)]
  data.frame(set = d$set[1], response = d$response[1], method = d$method[1], n = length(v),
             median = median(v), q10 = quantile(v, 0.1), q90 = quantile(v, 0.9),
             sign_top25 = median(d$sign_top25),
             corridoR_better = if (d$method[1] == "corridoR") NA else sum(ref > v),
             p = if (d$method[1] == "corridoR") NA else suppressWarnings(wilcox.test(ref, v, paired = TRUE, exact = FALSE)$p.value))
}))
write.csv(long, file.path(root, "summary.csv"), row.names = FALSE)
write.csv(bw, file.path(root, "bandwidths.csv"), row.names = FALSE)
saveRDS(x, file.path(root, "results.rds"))
print(long[order(long$set, long$response, -long$median), ], row.names = FALSE, digits = 3)
for (g in c("best_climate_bw", "best_wetness_bw", "best_terrain_bw", "best_cover_bw")) {
  cat("\n", g, "\n"); print(table(bw$set, bw[[g]]))
}
