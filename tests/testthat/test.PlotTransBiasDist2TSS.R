context("PlotTransBiasDist2TSS and PlotTransBiasDist2TSSToPDF")

test_that("PlotTransBiasDist2TSS function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh37.vcf")
  annotated.SBS.vcf <- AnnotateSBSVCF(list.of.vcfs$SBS.vcfs[[1]], "hg19", trans.ranges.GRCh37)
  out <- PlotTransBiasDist2TSS(annotated.SBS.vcf, "C>A")
  expect_equal(out, TRUE)

})

test_that("PlotTransBiasDist2TSSToPDF function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  list.of.vcfs <- ReadAndSplitStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh37.vcf")
  annotated.SBS.vcf <- AnnotateSBSVCF(list.of.vcfs$SBS.vcfs[[1]], "hg19", trans.ranges.GRCh37)
  out <- PlotTransBiasDist2TSSToPDF(annotated.SBS.vcf, "C>A",
                                    file = file.path(tempdir(), "test.pdf"))
  expect_equal(out, TRUE)
  
})