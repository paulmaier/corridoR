# A small GENEPOP file (6 meadows, 100 SNPs) cut from the example, for fast
# examples. Run after make_example.R, from the package root.
lines <- readLines("inst/extdata/example.gen")
pop_at <- which(lines == "Pop")
loci <- lines[2:(pop_at[1] - 1)]
keep_loci <- seq_len(100)
body <- unlist(lapply(c(3, 15, 28, 41, 55, 70), function(k) {
  end <- if (k < length(pop_at)) pop_at[k + 1] - 1 else length(lines)
  rows <- lines[(pop_at[k] + 1):end]
  g <- strsplit(sub("^[^,]*, *", "", rows), " ")
  c("Pop", paste0(sub(" *,.*$", "", rows), " , ", vapply(g, function(x) paste(x[keep_loci], collapse = " "), "")))
}))
writeLines(c("Small sample of the corridoR example", loci[keep_loci], body), "inst/extdata/example_small.gen")
