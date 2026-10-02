# Two populations, two loci, three diploid individuals; the expected allele
# counts are written out by hand below.
#   ind  pop  locus1  locus2
#   a1   A    1/2     1/1
#   a2   A    1/1     NA
#   b1   B    2/2     1/3
expect_counts <- function(pc) {
  expect_s3_class(pc, "popcounts")
  expect_equal(pc$pops, c("A", "B"))
  l1 <- pc$counts[[1]]; l2 <- pc$counts[[2]]
  expect_equal(unname(l1["A", ]), c(3L, 1L))
  expect_equal(unname(l1["B", ]), c(0L, 2L))
  expect_equal(unname(rowSums(l2)), c(2L, 2L))
  expect_equal(unname(l2["B", ]), c(1L, 1L))
}
tmp <- function(lines, ext) {
  f <- tempfile(fileext = ext)
  writeLines(lines, f)
  f
}

test_that("STRUCTURE, two rows per individual, with marker names", {
  f <- tmp(c("loc1 loc2",
             "a1 A 1 1", "a1 A 2 1",
             "a2 A 1 -9", "a2 A 1 -9",
             "b1 B 2 1", "b1 B 2 3"), ".str")
  pc <- read_genotypes(f)
  expect_counts(pc)
  expect_equal(names(pc$counts), c("loc1", "loc2"))
})

test_that("STRUCTURE, one row per individual, no header, extra columns", {
  f <- tmp(c("a1 A 0 1 2 1 1",
             "a2 A 0 1 1 -9 -9",
             "b1 B 1 2 2 1 3"), ".stru")
  expect_counts(read_genotypes(f, popflag = TRUE))
  # the same file without POPFLAG cannot pair its columns
  expect_error(read_genotypes(f), "cannot be paired")
})

test_that("STRUCTURE with populations supplied separately", {
  f <- tmp(c("a1 1 2 1 1", "a2 1 1 -9 -9", "b1 2 2 1 3"), ".str")
  pc <- read_genotypes(f, popdata = FALSE, onerowperind = TRUE,
                       pop = function(x) toupper(substr(x, 1, 1)))
  expect_equal(pc$pops, c("A", "B"))
  expect_equal(unname(pc$counts[[1]]["A", ]), c(3L, 1L))
})

test_that("GENEPOP with two- and three-digit codes and a blank title", {
  f2 <- tmp(c("", "loc1, loc2", "Pop", "a1 , 0102 0101", "a2 , 0101 0000",
              "Pop", "b1 , 0202 0103"), ".gen")
  pc <- read_genotypes(f2, pop_names = c("A", "B"))
  expect_counts(pc)
  f3 <- tmp(c("title", "loc1", "loc2", "POP", "a1, 001002 001001", "a2, 001001 000000",
              "pop", "b1, 002002 001003"), ".gen")
  expect_counts(read_genotypes(f3, pop_names = c("A", "B")))
  # by default a population is named after its first individual
  expect_equal(read_genotypes(f3)$pops, c("a1", "b1"))
})

test_that("VCF, phased and unphased, with missing calls and gzip", {
  v <- c("##fileformat=VCFv4.2",
         paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "a1", "a2", "b1", sep = "\t"),
         paste("1", "10", "loc1", "A", "G", ".", ".", ".", "GT:DP", "0/1:5", "0|0:7", "1/1:3", sep = "\t"),
         paste("1", "20", "loc2", "A", "G,T", ".", ".", ".", "GT", "0/0", "./.", "0/2", sep = "\t"))
  f <- tmp(v, ".vcf")
  pm <- data.frame(ind = c("a1", "a2", "b1"), pop = c("A", "A", "B"))
  expect_counts(read_genotypes(f, pop = pm))
  gz <- tempfile(fileext = ".vcf.gz")
  con <- gzfile(gz, "w"); writeLines(v, con); close(con)
  expect_counts(read_genotypes(gz, pop = c("A", "A", "B")))
  expect_error(read_genotypes(f), "Population assignments are needed")
})

test_that("GenAlEx", {
  f <- tmp(c("2,3,2,2,1", "Example,,A,B", "Ind,Pop,loc1,,loc2,",
             "a1,A,1,2,1,1", "a2,A,1,1,0,0", "b1,B,2,2,1,3"), ".csv")
  expect_counts(read_genotypes(f))
})

test_that("FSTAT", {
  f <- tmp(c("2 2 3 1", "loc1", "loc2", "1 12 11", "1 11 00", "2 22 13"), ".dat")
  pc <- read_genotypes(f, pop = c("A", "A", "B"))
  expect_counts(pc)
})

test_that("PLINK ped and raw", {
  ped <- tmp(c("A a1 0 0 0 -9 A G A A", "A a2 0 0 0 -9 A A 0 0", "B b1 0 0 0 -9 G G A C"), ".ped")
  writeLines(c("1 loc1 0 10", "1 loc2 0 20"), sub("\\.ped$", ".map", ped))
  expect_counts(read_genotypes(ped))
  raw <- tmp(c("FID IID PAT MAT SEX PHENOTYPE loc1_G loc2_C", "A a1 0 0 0 -9 1 0",
               "A a2 0 0 0 -9 0 NA", "B b1 0 0 0 -9 2 1"), ".raw")
  pc <- read_genotypes(raw)
  expect_equal(unname(pc$counts[[1]]["A", ]), c(3L, 1L))
  expect_equal(unname(pc$counts[[2]]["B", ]), c(1L, 1L))
})

test_that("the example files agree across formats", {
  ex <- corridor_example()
  a <- read_genotypes(ex$files$structure)
  b <- read_genotypes(ex$files$genepop, pop = function(x) sub("_.*", "", x))
  v <- read_genotypes(ex$files$vcf, pop = ex$files$popmap)
  expect_equal(pairwise_fst(a), pairwise_fst(b))
  expect_equal(pairwise_fst(a), pairwise_fst(v))
})

test_that("adegenet objects convert", {
  skip_if_not_installed("adegenet")
  ex <- corridor_example()
  gi <- suppressWarnings(adegenet::read.genepop(ex$files$genepop, quiet = TRUE))
  pc <- as_popcounts(gi)
  ref <- read_genotypes(ex$files$genepop)
  expect_equal(length(pc$counts), length(ref$counts))
  expect_equal(unname(pairwise_fst(pc)), unname(pairwise_fst(ref)))
})
