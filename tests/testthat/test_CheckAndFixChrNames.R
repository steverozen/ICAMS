context("CheckAndFixChrNames")

test_that("CheckAndFixChrNames for human genome 37", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  
  make.input <- function(x) {
    return(data.frame(CHROM = x, stringsAsFactors = FALSE))
  }
  
  in.vcf <- make.input(c("1", "2")) 
  
  expect_equal(
    CheckAndFixChrNames(
      vcf.df = in.vcf,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5),
    c("1", "2"))
  
  in.vcf1 <- make.input(c("chr1", "chr2")) 
  expect_equal(
    CheckAndFixChrNames(
      vcf.df = in.vcf1,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5),
    c("1", "2"))
  
  in.vcf2 <- make.input(c("1", "23"))
  expect_equal(
    CheckAndFixChrNames(
      vcf.df = in.vcf2,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5),
    c("1", "X"))
  
  in.vcf3 <- make.input(c("chr1", "chr23"))
  expect_equal(
    CheckAndFixChrNames(
      vcf.df = in.vcf3,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5),
    c("1", "X"))
  
  in.vcf4 <- make.input(c("23", "24"))
  expect_equal(
    CheckAndFixChrNames(
      vcf.df = in.vcf4,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5),
    c("X", "Y"))
  
  in.vcf5 <- make.input(c("chr23", "chr24"))
  expect_equal(
    CheckAndFixChrNames(
      vcf.df = in.vcf5,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5),
    c("X", "Y"))
  
  in.vcf6 <- make.input(c("23", "X"))
  expect_error(
    CheckAndFixChrNames(
      vcf.df = in.vcf6,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5))
  
  in.vcf7 <- make.input(c("chr23", "chrX"))
  expect_error(
    CheckAndFixChrNames(
      vcf.df = in.vcf7,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5))
  
  in.vcf8 <- make.input(c("24", "Y"))
  expect_error(
    CheckAndFixChrNames(
      vcf.df = in.vcf8,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5))
  
  in.vcf9 <- make.input(c("chr24", "chrY"))
  expect_error(
    CheckAndFixChrNames(
      vcf.df = in.vcf9,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5))
  
  bad.vcf <- make.input(c("1", "chr2"))
  expect_error(
    CheckAndFixChrNames(
      vcf.df = bad.vcf,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5))
}
)

test_that("CheckAndFixChrNames for human genome 38", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.UCSC.hg38"))
  
  make.input <- function(x) {
    return(data.frame(CHROM = x, stringsAsFactors = FALSE))
  }
  
  in.vcf <- make.input(c("1", "2")) 
  expect_equal(
    CheckAndFixChrNames(
      vcf.df = in.vcf,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38),
    c("chr1", "chr2"))
  
  in.vcf1 <- make.input(c("chr1", "chr2")) 
  expect_equal(
    CheckAndFixChrNames(
      vcf.df = in.vcf1,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38),
    c("chr1", "chr2"))
  
  in.vcf2 <- make.input(c("1", "23"))
  expect_equal(
    CheckAndFixChrNames(
      vcf.df = in.vcf2,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38),
    c("chr1", "chrX"))
  
  in.vcf3 <- make.input(c("chr1", "chr23"))
  expect_equal(
    CheckAndFixChrNames(
      vcf.df = in.vcf3,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38),
    c("chr1", "chrX"))
  
  in.vcf4 <- make.input(c("23", "24"))
  expect_equal(
    CheckAndFixChrNames(
      vcf.df = in.vcf4,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38),
    c("chrX", "chrY"))
  
  in.vcf5 <- make.input(c("chr23", "chr24"))
  expect_equal(
    CheckAndFixChrNames(
      vcf.df = in.vcf4,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38),
    c("chrX", "chrY"))
  
  in.vcf6 <- make.input(c("23", "X"))
  expect_error(
    CheckAndFixChrNames(
      vcf.df = in.vcf6,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))
  
  in.vcf7 <- make.input(c("chr23", "chrX"))
  expect_error(
    CheckAndFixChrNames(
      vcf.df = in.vcf7,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))
  
  in.vcf8 <- make.input(c("24", "Y"))
  expect_error(
    CheckAndFixChrNames(
      vcf.df = in.vcf8,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))
  
  in.vcf9 <- make.input(c("chr24", "chrY"))
  expect_error(
    CheckAndFixChrNames(
      vcf.df = in.vcf9,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))
  
  bad.vcf <- make.input(c("1", "chr2"))
  expect_error(
    CheckAndFixChrNames(
      vcf.df = bad.vcf,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))
}
)

test_that("CheckAndFixChrNames for mouse genome", {
  skip_if("" == system.file(package = "BSgenome.Mmusculus.UCSC.mm10"))
  stopifnot(requireNamespace("BSgenome.Mmusculus.UCSC.mm10"))
  
  make.input <- function(x) {
    return(data.frame(CHROM = x, stringsAsFactors = FALSE))
  }
  
  in.vcf <- make.input(c("1", "2")) 
  expect_equal(
    CheckAndFixChrNames(
      vcf.df = in.vcf,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10),
    c("chr1", "chr2"))
  
  in.vcf1 <- make.input(c("chr1", "chr2")) 
  expect_equal(
    CheckAndFixChrNames(
      vcf.df = in.vcf1,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10),
    c("chr1", "chr2"))
  
  in.vcf2 <- make.input(c("1", "20"))
  expect_equal(
    CheckAndFixChrNames(
      vcf.df = in.vcf2,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10),
    c("chr1", "chrX"))
  
  in.vcf3 <- make.input(c("chr1", "chr20"))
  expect_equal(
    CheckAndFixChrNames(
      vcf.df = in.vcf3,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10),
    c("chr1", "chrX"))
  
  in.vcf4 <- make.input(c("20", "21"))
  expect_equal(
    CheckAndFixChrNames(
      vcf.df = in.vcf4,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10),
    c("chrX", "chrY"))
  
  in.vcf5 <- make.input(c("chr20", "chr21"))
  expect_equal(
    CheckAndFixChrNames(
      vcf.df = in.vcf4,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10),
    c("chrX", "chrY"))
  
  in.vcf6 <- make.input(c("20", "X"))
  expect_error(
    CheckAndFixChrNames(
      vcf.df = in.vcf6,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10))
  
  in.vcf7 <- make.input(c("chr20", "chrX"))
  expect_error(
    CheckAndFixChrNames(
      vcf.df = in.vcf7,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10))
  
  in.vcf8 <- make.input(c("21", "Y"))
  expect_error(
    CheckAndFixChrNames(
      vcf.df = in.vcf8,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10))
  
  in.vcf9 <- make.input(c("chr21", "chrY"))
  expect_error(
    CheckAndFixChrNames(
      vcf.df = in.vcf9,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10))
  
  bad.vcf <- make.input(c("1", "chr2"))
  expect_error(
    CheckAndFixChrNames(
      vcf.df = bad.vcf,
      ref.genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10))
}
)

