context("CheckAndFixChrNamesForTransRanges")

test_that("CheckAndFixChrNamesForTransRanges for human genome 37", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  
  make.input <- function(x) {
    return(data.frame(CHROM = x, stringsAsFactors = FALSE))
  }
  
  in.vcf <- make.input(c("1", "2")) 
  
  trans.ranges <- trans.ranges.GRCh37
  
  expect_equal(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5),
    trans.ranges$chrom)
  
  in.vcf1 <- make.input(c("chr1", "chr2")) 
  expect_equal(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf1,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5),
    paste0("chr", trans.ranges$chrom))
  
  in.vcf2 <- make.input(c("1", "23"))
  retval.to.check2 <- as.character(trans.ranges$chrom)
  retval.to.check2[retval.to.check2 == "X"] <- "23"
  expect_equal(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf2,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5),
    retval.to.check2)
  
  in.vcf3 <- make.input(c("chr1", "chr23"))
  retval.to.check3 <- paste0("chr", as.character(trans.ranges$chrom))
  retval.to.check3[retval.to.check3 == "chrX"] <- "chr23"
  expect_equal(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf3,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5),
    retval.to.check3)
  
  in.vcf4 <- make.input(c("23", "24"))
  retval.to.check4 <- as.character(trans.ranges$chrom)
  retval.to.check4[retval.to.check4 == "X"] <- "23"
  retval.to.check4[retval.to.check4 == "Y"] <- "24"
  expect_equal(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf4,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5),
    retval.to.check4)
  
  in.vcf5 <- make.input(c("chr23", "chr24"))
  retval.to.check5 <- paste0("chr", as.character(trans.ranges$chrom))
  retval.to.check5[retval.to.check5 == "chrX"] <- "chr23"
  retval.to.check5[retval.to.check5 == "chrY"] <- "chr24"
  expect_equal(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf5,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5),
    retval.to.check5)
  
  in.vcf6 <- make.input(c("23", "X"))
  expect_error(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf6,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5))
  
  in.vcf7 <- make.input(c("chr23", "chrX"))
  expect_error(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf7,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5))
  
  in.vcf8 <- make.input(c("24", "Y"))
  expect_error(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf8,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5))
  
  in.vcf9 <- make.input(c("chr24", "chrY"))
  expect_error(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf9,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5))
  
  bad.vcf <- make.input(c("1", "chr2"))
  expect_error(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = bad.vcf,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5))
}
)

test_that("CheckAndFixChrNamesForTransRanges for human genome 38", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.UCSC.hg38"))
  
  make.input <- function(x) {
    return(data.frame(CHROM = x, stringsAsFactors = FALSE))
  }
  
  in.vcf <- make.input(c("1", "2")) 
  
  trans.ranges <- trans.ranges.GRCh38
  
  expect_equal(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38),
    trans.ranges$chrom)
  
  in.vcf1 <- make.input(c("chr1", "chr2")) 
  expect_equal(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf1,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38),
    paste0("chr", trans.ranges$chrom))
  
  in.vcf2 <- make.input(c("1", "23"))
  retval.to.check2 <- as.character(trans.ranges$chrom)
  retval.to.check2[retval.to.check2 == "X"] <- "23"
  expect_equal(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf2,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38),
    retval.to.check2)
  
  in.vcf3 <- make.input(c("chr1", "chr23"))
  retval.to.check3 <- paste0("chr", as.character(trans.ranges$chrom))
  retval.to.check3[retval.to.check3 == "chrX"] <- "chr23"
  expect_equal(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf3,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38),
    retval.to.check3)
  
  in.vcf4 <- make.input(c("23", "24"))
  retval.to.check4 <- as.character(trans.ranges$chrom)
  retval.to.check4[retval.to.check4 == "X"] <- "23"
  retval.to.check4[retval.to.check4 == "Y"] <- "24"
  expect_equal(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf4,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38),
    retval.to.check4)
  
  in.vcf5 <- make.input(c("chr23", "chr24"))
  retval.to.check5 <- paste0("chr", as.character(trans.ranges$chrom))
  retval.to.check5[retval.to.check5 == "chrX"] <- "chr23"
  retval.to.check5[retval.to.check5 == "chrY"] <- "chr24"
  expect_equal(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf5,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38),
    retval.to.check5)
  
  in.vcf6 <- make.input(c("23", "X"))
  expect_error(
    CheckAndFixChrNamesForTransRanges(
      vcf.df = in.vcf6,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))
  
  in.vcf7 <- make.input(c("chr23", "chrX"))
  expect_error(
    CheckAndFixChrNamesForTransRanges(
      vcf.df = in.vcf7,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))
  
  in.vcf8 <- make.input(c("24", "Y"))
  expect_error(
    CheckAndFixChrNamesForTransRanges(
      vcf.df = in.vcf8,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))
  
  in.vcf9 <- make.input(c("chr24", "chrY"))
  expect_error(
    CheckAndFixChrNamesForTransRanges(
      vcf.df = in.vcf9,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))
  
  bad.vcf <- make.input(c("1", "chr2"))
  expect_error(
    CheckAndFixChrNamesForTransRanges(
      vcf.df = bad.vcf,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))
}
)

test_that("CheckAndFixChrNamesForTransRanges for mouse genome", {
  skip_if("" == system.file(package = "BSgenome.Mmusculus.UCSC.mm10"))
  stopifnot(requireNamespace("BSgenome.Mmusculus.UCSC.mm10"))
  
  make.input <- function(x) {
    return(data.frame(CHROM = x, stringsAsFactors = FALSE))
  }
  
  in.vcf <- make.input(c("1", "2")) 
  
  trans.ranges <- trans.ranges.GRCm38
  
  expect_equal(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10),
    trans.ranges$chrom)
  
  in.vcf1 <- make.input(c("chr1", "chr2")) 
  expect_equal(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf1,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10),
    paste0("chr", trans.ranges$chrom))
  
  in.vcf2 <- make.input(c("1", "20"))
  retval.to.check2 <- as.character(trans.ranges$chrom)
  retval.to.check2[retval.to.check2 == "X"] <- "20"
  expect_equal(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf2,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10),
    retval.to.check2)
  
  in.vcf3 <- make.input(c("chr1", "chr20"))
  retval.to.check3 <- paste0("chr", as.character(trans.ranges$chrom))
  retval.to.check3[retval.to.check3 == "chrX"] <- "chr20"
  expect_equal(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf3,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10),
    retval.to.check3)
  
  in.vcf4 <- make.input(c("20", "21"))
  retval.to.check4 <- as.character(trans.ranges$chrom)
  retval.to.check4[retval.to.check4 == "X"] <- "20"
  retval.to.check4[retval.to.check4 == "Y"] <- "21"
  expect_equal(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf4,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10),
    retval.to.check4)
  
  in.vcf5 <- make.input(c("chr20", "chr21"))
  retval.to.check5 <- paste0("chr", as.character(trans.ranges$chrom))
  retval.to.check5[retval.to.check5 == "chrX"] <- "chr20"
  retval.to.check5[retval.to.check5 == "chrY"] <- "chr21"
  expect_equal(
    CheckAndFixChrNamesForTransRanges(
      trans.ranges = trans.ranges,
      vcf.df = in.vcf5,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10),
    retval.to.check5)
  
  in.vcf6 <- make.input(c("20", "X"))
  expect_error(
    CheckAndFixChrNamesForTransRanges(
      vcf.df = in.vcf6,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10))
  
  in.vcf7 <- make.input(c("chr20", "chrX"))
  expect_error(
    CheckAndFixChrNamesForTransRanges(
      vcf.df = in.vcf7,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10))
  
  in.vcf8 <- make.input(c("21", "Y"))
  expect_error(
    CheckAndFixChrNamesForTransRanges(
      vcf.df = in.vcf8,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10))
  
  in.vcf9 <- make.input(c("chr21", "chrY"))
  expect_error(
    CheckAndFixChrNamesForTransRanges(
      vcf.df = in.vcf9,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10))
  
  bad.vcf <- make.input(c("1", "chr2"))
  expect_error(
    CheckAndFixChrNamesForTransRanges(
      vcf.df = bad.vcf,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10))
}
)

