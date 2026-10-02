# Study-area figure for the tutorial: meadows on the terrain, with a lineage
# tree as the legend. Run from the package root.
devtools::load_all(quiet = TRUE); library(terra); library(sf)
ex <- corridor_example()
lc <- c(North = "#ff3939", East = "#0070ff", West = "#000000", South = "#4ce600")
hill <- shade(terrain(ex$dem * 3, "slope", unit = "radians"), terrain(ex$dem, "aspect", unit = "radians"), 35, 315)
dem_cols <- colorRampPalette(c("#5b7f4a", "#9fae6a", "#d8cf98", "#c9b28a", "#ece6dc", "#ffffff"))(80)

png("man/figures/study_area.png", width = 2200, height = 1250, res = 200)
layout(matrix(1:2, 1), widths = c(4, 1.25))
par(mar = c(1, 1, 2, 1))
plot(hill, col = grey(0:100 / 100), legend = FALSE, axes = FALSE, main = "Simulated range: meadows and lineages")
plot(ex$dem, col = dem_cols, alpha = 0.8, add = TRUE,
     plg = list(title = "Elevation (m)", x = "bottomleft", horiz = TRUE, cex = 0.7))
plot(st_geometry(ex$sites), add = TRUE, pch = 21, bg = lc[ex$sites$lineage], col = "white", cex = 1.6, lwd = 1.1)
text(c(4000, 4000), c(30500, 13800), c("north canyon", "south canyon"), col = "grey15", cex = 0.7, adj = 0, font = 3)
text(50500, 39000, "crest", col = "grey15", cex = 0.75, font = 3)

# lineage tree from lineage_tmrca: (North, ((East, West), South))
tm <- ex$lineage_tmrca
age <- function(a, b) tm$tmrca[(tm$lineage1 == a & tm$lineage2 == b) | (tm$lineage1 == b & tm$lineage2 == a)]
root <- age("North", "West"); ew <- age("East", "West"); s_node <- age("West", "South")
yy <- c(North = 4, East = 3, West = 2, South = 1)
par(mar = c(4, 0.5, 2, 0.5))
plot(NA, xlim = c(root * 1.12, -root * 0.55), ylim = c(0.5, 4.6), axes = FALSE, xlab = "", ylab = "",
     main = "Lineages", cex.main = 1)
seg <- function(x0, y0, x1, y1, ...) segments(x0, y0, x1, y1, lwd = 3, lend = 1, ...)
y_ew <- mean(yy[c("East", "West")]); y_s <- mean(c(y_ew, yy["South"])); y_root <- mean(c(yy["North"], y_s))
seg(ew, yy["East"], ew, yy["West"]); seg(s_node, y_ew, s_node, yy["South"]); seg(root, y_s, root, yy["North"])
seg(root, y_s, s_node, y_s); seg(s_node, y_ew, ew, y_ew)
for (l in names(yy)) {
  from <- switch(l, North = root, South = s_node, ew)
  seg(from, yy[l], 0, yy[l], col = lc[l])
  points(0, yy[l], pch = 21, bg = lc[l], col = "white", cex = 2)
  text(-root * 0.06, yy[l], l, adj = 0, cex = 0.9)
}
axis(1, at = seq(0, root, by = 200), labels = seq(0, root, by = 200), cex.axis = 0.7)
mtext("Thousand years ago", 1, line = 2.3, cex = 0.7)
dev.off()
