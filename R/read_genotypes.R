#' Read genotypes from a file
#'
#' Reads the common population genetic file formats into a [popcounts][as_popcounts]
#' object (allele counts per population and locus), which is what
#' [pairwise_fst()] and [directional_gst()] use.
#'
#' @section Supported formats:
#' \describe{
#' \item{STRUCTURE (`.str`, `.stru`, `.structure`)}{Both layouts: one row per
#'   individual (two columns per locus) or two rows per individual (one column
#'   per locus). Optional leading columns, in STRUCTURE's order: a label, the
#'   population (POPDATA), POPFLAG, LOCDATA, PHENOTYPE and any extra columns.
#'   An optional first row of marker names is detected automatically. Alleles
#'   may be any integer codes; missing data are `-9` by default. When the file
#'   has a marker-name row, the number of leading columns is worked out from
#'   it, otherwise set the arguments below.}
#' \item{GENEPOP (`.gen`, `.genepop`)}{Title line, loci one per line or
#'   comma-separated, `Pop` separators, and `name , genotypes` rows. Two- or
#'   three-digit allele codes (detected from the genotype width), haploid or
#'   diploid; `0`, `00`, `000` and so on are missing.}
#' \item{VCF (`.vcf`, `.vcf.gz`)}{Any ploidy, phased or unphased, biallelic or
#'   multiallelic. VCF has no populations, so supply `pop`.}
#' \item{GenAlEx (`.csv`)}{The standard codominant layout: a first row with the
#'   number of loci, individuals and populations, a second row with population
#'   names, a header row, then sample, population and two columns per locus.}
#' \item{FSTAT (`.dat`)}{The standard header (populations, loci, maximum allele,
#'   digits per allele), locus names, then population and genotype columns.}
#' \item{PLINK (`.ped` with `.map`, or `.raw` from `--recode A`)}{The family ID
#'   is used as the population unless `pop` is given. Convert binary `.bed`
#'   files with `plink --recode A` or to VCF first.}
#' }
#' Objects already in R (\pkg{adegenet} `genind`, `genpop`, `genlight`, or a
#' \pkg{vcfR} object) go through [as_popcounts()] instead.
#'
#' @param file Path to the genotype file.
#' @param format One of `"auto"` (guess from the extension), `"structure"`,
#'   `"genepop"`, `"vcf"`, `"genalex"`, `"fstat"`, `"plink"`.
#' @param pop Optional population assignment that overrides the one in the
#'   file: a vector with one value per individual in file order, a
#'   two-column data.frame or whitespace-separated file (individual,
#'   population), or a function that takes individual names and returns
#'   populations (for example `function(x) sub("_.*", "", x)`).
#' @param missing STRUCTURE only: codes for missing alleles.
#' @param onerowperind,label,popdata,popflag,locdata,phenotype,extracols
#'   STRUCTURE only. `NULL` means detect: the row layout is detected from
#'   repeated labels, and the leading columns from the marker-name row when
#'   there is one. Otherwise the defaults are a label and a population column,
#'   which is the most common layout.
#' @param pop_names GENEPOP only: names for the populations, in file order.
#'   By default each population takes the name of its first individual, which
#'   is the GENEPOP convention.
#' @return A `popcounts` object.
#' @examples
#' ex <- corridor_example()
#' read_genotypes(ex$files$structure)
#' read_genotypes(ex$files$genepop)
#' read_genotypes(ex$files$vcf, pop = ex$files$popmap)
#' @export
read_genotypes <- function(file, format = c("auto", "structure", "genepop", "vcf", "genalex",
                                            "fstat", "plink"),
                           pop = NULL, missing = c("-9"), onerowperind = NULL, label = NULL,
                           popdata = NULL, popflag = FALSE, locdata = FALSE, phenotype = FALSE,
                           extracols = 0, pop_names = NULL) {
  if (!file.exists(file)) stop("Cannot find ", file, call. = FALSE)
  format <- match.arg(format)
  if (format == "auto") format <- guess_format(file)
  switch(format,
         structure = read_structure(file, pop, missing, onerowperind, label, popdata,
                                    popflag, locdata, phenotype, extracols),
         genepop = read_genepop(file, pop, pop_names),
         vcf = read_vcf(file, pop),
         genalex = read_genalex(file, pop),
         fstat = read_fstat(file, pop),
         plink = read_plink(file, pop))
}

guess_format <- function(file) {
  f <- tolower(basename(file))
  f <- sub("\\.gz$", "", f)
  ext <- sub(".*\\.", "", f)
  out <- switch(ext, str = , stru = , structure = "structure", gen = , genepop = "genepop",
                vcf = "vcf", csv = "genalex", dat = "fstat", ped = , raw = "plink", NA)
  if (is.na(out)) {
    stop("Cannot tell the format of ", basename(file), " from its extension; set `format`.",
         call. = FALSE)
  }
  out
}

# ---- STRUCTURE ------------------------------------------------------------------

read_structure <- function(file, pop, missing, onerowperind, label, popdata, popflag, locdata,
                           phenotype, extracols) {
  lines <- trimws(readLines(file, warn = FALSE))
  lines <- lines[nzchar(lines)]
  tok <- strsplit(lines, "[[:space:],]+")
  len <- lengths(tok)
  body_len <- as.integer(names(which.max(table(len))))
  header <- NULL
  if (len[1] < body_len) {                       # marker names (and maybe map distances)
    header <- tok[[1]]
    drop <- 1L
    if (length(len) > 1 && len[2] < body_len) drop <- 2L
    tok <- tok[-seq_len(drop)]; len <- len[-seq_len(drop)]
  }
  if (any(len != body_len)) {
    bad <- which(len != body_len)[1]
    stop("Row ", bad, " of the genotype block has ", len[bad], " fields; the others have ",
         body_len, ". Check for spaces inside sample names.", call. = FALSE)
  }
  m <- do.call(rbind, tok)

  if (is.null(label)) label <- !all(grepl("^-?[0-9]+$", m[, 1]))
  if (is.null(onerowperind)) {
    onerowperind <- !(label && nrow(m) %% 2 == 0 &&
                        all(m[seq(1, nrow(m), 2), 1] == m[seq(2, nrow(m), 2), 1]))
  }
  if (is.null(popdata)) {
    if (!is.null(header)) {
      n_meta <- body_len - length(header) * if (onerowperind) 2L else 1L
      popdata <- n_meta - label - popflag - locdata - phenotype - extracols >= 1
    } else {
      popdata <- TRUE
    }
  }
  n_meta <- label + popdata + popflag + locdata + phenotype + extracols
  geno <- m[, -seq_len(n_meta), drop = FALSE]
  nloc <- if (onerowperind) ncol(geno) / 2 else ncol(geno)
  if (nloc != round(nloc)) {
    stop("With ", n_meta, " leading columns and one row per individual, the genotype columns ",
         "cannot be paired into loci. Check `label`, `popdata`, `popflag`, `locdata`, ",
         "`phenotype` and `extracols`.", call. = FALSE)
  }
  if (!is.null(header) && length(header) != nloc) {
    stop("The marker-name row lists ", length(header), " loci, but the data have ", nloc,
         " given the leading columns. Check the STRUCTURE options.", call. = FALSE)
  }
  geno[geno %in% missing] <- NA
  if (onerowperind) {
    ind <- if (label) m[, 1] else paste0("ind", seq_len(nrow(m)))
    arr <- array(NA_character_, c(nrow(m), nloc, 2))
    arr[, , 1] <- geno[, seq(1, ncol(geno), 2)]
    arr[, , 2] <- geno[, seq(2, ncol(geno), 2)]
    meta_rows <- seq_len(nrow(m))
  } else {
    r1 <- seq(1, nrow(m), 2)
    ind <- if (label) m[r1, 1] else paste0("ind", seq_along(r1))
    arr <- array(NA_character_, c(length(r1), nloc, 2))
    arr[, , 1] <- geno[r1, ]
    arr[, , 2] <- geno[r1 + 1, ]
    meta_rows <- r1
  }
  if (is.null(pop)) {
    if (!popdata) stop("This STRUCTURE file has no population column; supply `pop`.", call. = FALSE)
    pop <- m[meta_rows, label + 1]
  }
  loci <- if (is.null(header)) paste0("locus", seq_len(nloc)) else header
  calls_to_popcounts(arr, resolve_pop(pop, ind), loci)
}

# ---- GENEPOP ------------------------------------------------------------------

read_genepop <- function(file, pop, pop_names) {
  lines <- trimws(readLines(file, warn = FALSE))
  lines <- c(lines[1], lines[-1][nzchar(lines[-1])])   # the title line may be blank
  is_pop <- grepl("^pop$", lines, ignore.case = TRUE) | grepl("^pop[[:space:],]", lines, ignore.case = TRUE)
  first <- which(is_pop)[1]
  if (is.na(first)) stop("No `Pop` line found; is this a GENEPOP file?", call. = FALSE)
  loci <- unlist(strsplit(lines[2:(first - 1)], "[[:space:]]*,[[:space:]]*"))
  loci <- trimws(loci[nzchar(trimws(loci))])
  body <- lines[first:length(lines)]
  block <- cumsum(grepl("^pop", body, ignore.case = TRUE) & !grepl(",", body))
  rows <- body[!(grepl("^pop", body, ignore.case = TRUE) & !grepl(",", body))]
  grp <- block[!(grepl("^pop", body, ignore.case = TRUE) & !grepl(",", body))]
  ind <- trimws(sub(",.*$", "", rows))
  gts <- strsplit(trimws(sub("^[^,]*,", "", rows)), "[[:space:]]+")
  if (any(lengths(gts) != length(loci))) {
    stop("Some individuals do not have ", length(loci), " genotypes.", call. = FALSE)
  }
  g <- do.call(rbind, gts)
  width <- max(nchar(g))
  digits <- if (width %in% c(4, 6)) width / 2 else width        # diploid 2/3-digit, else haploid
  ploidy <- width / digits
  arr <- array(NA_character_, c(nrow(g), ncol(g), ploidy))
  for (k in seq_len(ploidy)) {
    a <- substr(g, (k - 1) * digits + 1, k * digits)
    a[!grepl("[1-9]", a)] <- NA
    arr[, , k] <- a
  }
  if (is.null(pop)) {
    nm <- if (is.null(pop_names)) ind[!duplicated(grp)] else pop_names
    if (length(nm) != max(grp)) stop("`pop_names` needs ", max(grp), " names.", call. = FALSE)
    pop <- nm[grp]
  }
  calls_to_popcounts(arr, resolve_pop(pop, ind), loci)
}

# ---- VCF ------------------------------------------------------------------

read_vcf <- function(file, pop) {
  con <- if (grepl("\\.gz$", file)) gzfile(file) else file(file)
  lines <- readLines(con, warn = FALSE)
  close(con)
  lines <- lines[!startsWith(lines, "##")]
  hdr <- strsplit(sub("^#", "", lines[1]), "\t")[[1]]
  f <- strsplit(lines[-1], "\t", fixed = TRUE)
  f <- do.call(rbind, f)
  samples <- hdr[-(1:9)]
  fmt <- strsplit(f[, 9], ":", fixed = TRUE)
  gt_i <- vapply(fmt, function(z) match("GT", z), integer(1))
  calls <- f[, -(1:9), drop = FALSE]
  if (!all(gt_i == 1L)) {
    for (i in which(gt_i != 1L)) calls[i, ] <- vapply(strsplit(calls[i, ], ":", fixed = TRUE), `[`, "", gt_i[i])
  } else {
    calls[] <- sub(":.*$", "", calls)
  }
  ids <- f[, 3]
  ids[ids == "."] <- paste(f[ids == ".", 1], f[ids == ".", 2], sep = "_")
  colnames(calls) <- samples
  vcf_gt_counts(calls, ids, resolve_pop(pop, samples))
}

# ---- GenAlEx ---------------------------------------------------------------

read_genalex <- function(file, pop) {
  x <- utils::read.csv(file, header = FALSE, colClasses = "character", strip.white = TRUE)
  nloc <- as.integer(x[1, 1]); nind <- as.integer(x[1, 2])
  if (is.na(nloc) || is.na(nind)) stop("The first row should hold the number of loci and individuals.", call. = FALSE)
  loci <- as.character(x[3, seq(3, by = 2, length.out = nloc)])
  d <- x[4:(3 + nind), , drop = FALSE]
  g <- as.matrix(d[, 3:(2 + 2 * nloc)])
  g[g %in% c("0", "-9", "", "?")] <- NA
  arr <- array(NA_character_, c(nind, nloc, 2))
  arr[, , 1] <- g[, seq(1, 2 * nloc, 2)]
  arr[, , 2] <- g[, seq(2, 2 * nloc, 2)]
  if (is.null(pop)) pop <- d[[2]]
  calls_to_popcounts(arr, resolve_pop(pop, d[[1]]), loci)
}

# ---- FSTAT ----------------------------------------------------------------------

read_fstat <- function(file, pop) {
  lines <- trimws(readLines(file, warn = FALSE))
  lines <- lines[nzchar(lines)]
  h <- as.integer(strsplit(lines[1], "[[:space:]]+")[[1]])
  nloc <- h[2]; digits <- h[4]
  loci <- lines[2:(1 + nloc)]
  g <- do.call(rbind, strsplit(lines[-(1:(1 + nloc))], "[[:space:]]+"))
  pops <- g[, 1]; g <- g[, -1, drop = FALSE]
  g <- formatC(as.integer(g), width = 2 * digits, flag = "0")
  g <- matrix(g, nrow = length(pops))
  arr <- array(NA_character_, c(nrow(g), nloc, 2))
  for (k in 1:2) {
    a <- substr(g, (k - 1) * digits + 1, k * digits)
    a[!grepl("[1-9]", a)] <- NA
    arr[, , k] <- a
  }
  if (is.null(pop)) pop <- pops
  calls_to_popcounts(arr, resolve_pop(pop, paste0("ind", seq_along(pops))), loci)
}

# ---- PLINK ----------------------------------------------------------------------

read_plink <- function(file, pop) {
  if (grepl("\\.raw$", file)) {
    x <- utils::read.table(file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    g <- as.matrix(x[, -(1:6), drop = FALSE])
    colnames(g) <- sub("_[^_]*$", "", colnames(g))
    if (is.null(pop)) pop <- x$FID
    return(dosage_counts(g, resolve_pop(pop, x$IID), 2L))
  }
  x <- utils::read.table(file, header = FALSE, colClasses = "character")
  map <- sub("\\.ped$", ".map", file)
  loci <- if (file.exists(map)) utils::read.table(map, colClasses = "character")[[2]] else
    paste0("snp", seq_len((ncol(x) - 6) / 2))
  g <- as.matrix(x[, -(1:6), drop = FALSE])
  g[g %in% c("0", "-9", "N")] <- NA
  nloc <- ncol(g) / 2
  arr <- array(NA_character_, c(nrow(g), nloc, 2))
  arr[, , 1] <- g[, seq(1, ncol(g), 2)]
  arr[, , 2] <- g[, seq(2, ncol(g), 2)]
  if (is.null(pop)) pop <- x[[1]]
  calls_to_popcounts(arr, resolve_pop(pop, x[[2]]), loci)
}
