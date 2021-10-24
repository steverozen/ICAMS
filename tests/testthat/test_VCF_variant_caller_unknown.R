context("Test processing VCF with unknown variant caller")

test_that("Test processing VCF with unknown variant caller", {
  #skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  #stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  file <- "testdata/SBS.GRCh37.variantcaller.unknown.vcf"
  
  # Do not merge adjacent SBSs into DBS
  list.of.vcfs1 <- ReadAndSplitVCFs(file, filter.status = "PASS")
  
  # Merge adjacent SBSs into DBS
  list.of.vcfs2 <- ReadAndSplitVCFs(file, get.vaf.function = function(x){
    x$VAF <- 0.5
    x$read.depth <- NA
    return(x)
  }, filter.status = "PASS")
  expect_equal(nrow(list.of.vcfs1$DBS[[1]]), 0)
  expect_equal(nrow(list.of.vcfs2$DBS[[1]]), 18)
})