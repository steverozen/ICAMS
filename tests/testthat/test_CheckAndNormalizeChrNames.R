context("CheckAndNormalizeChrNames")

test_that("All CheckAndNormalizeChrNames", {
  
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.UCSC.hg38"))
  
  make.input <- function(x) {
    return(data.frame(CHROM = x, stringsAsFactors = FALSE))
  }
  
  in.vcf <- make.input(c("1", "2")) 
  
  expect_equal(
    CheckAndNormalizeChrNames(
      vcf.df = in.vcf,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5),
    c("1", "2"))
  
  expect_equal(
    CheckAndNormalizeChrNames(
      vcf.df = in.vcf,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38),
    c("chr1", "chr2"))
  
  in.vcf <- make.input(c("1", "23"))
  expect_equal(
    CheckAndNormalizeChrNames(
      vcf.df = in.vcf,
      ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5),
    c("1", "X"))
  
  expect_equal(
    CheckAndNormalizeChrNames(
      vcf.df = in.vcf,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38),
    c("chr1", "chrX"))
  
  bad.vcf <- make.input(c("1", "chr2"))
  expect_error(
    CheckAndNormalizeChrNames(
      vcf.df = bad.vcf,
      ref.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))

  
}
)