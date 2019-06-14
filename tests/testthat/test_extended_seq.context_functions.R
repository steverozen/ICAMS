context("extended seq.context functions")

test_that("extended seq.context functions for Mutect GRCh38 vcf", {
  vcf <- ReadAndSplitMutectVCFs("testdata/Mutect.GRCh38.vcf")
  sbs.vcf <- vcf$SBS[[1]]
  mat <- CreateOnePPMFromSBSVCF(sbs.vcf, ref.genome = "GRCh38",
                                seq.context.width = 10)
  out <- PlotPPM(mat, title = "ExtendedSeqContext_21bases")
  out1 <- 
    PlotPPMToPdf(list(mat), titles = "ExtendedSeqContext_21bases",
                 file = paste0(tempdir(), "\\Extended.seq.context.21bases.pdf"))
  expect_equal(out, TRUE)
  expect_equal(out1, TRUE)
  unlink(paste0(tempdir(), "\\Extended.seq.context.21bases.pdf"))
})

test_that("extended seq.context functions for Mutect GRCh37 vcf", {
  vcf <- ReadAndSplitMutectVCFs("testdata/Mutect.GRCh37.vcf")
  sbs.vcf <- vcf$SBS[[1]]
  mat <- CreateOnePPMFromSBSVCF(sbs.vcf, ref.genome = "GRCh37",
                                seq.context.width = 10)
  out <- PlotPPM(mat, title = "ExtendedSeqContext_21bases")
  out1 <- 
    PlotPPMToPdf(list(mat), titles = "ExtendedSeqContext_21bases",
                 file = paste0(tempdir(), "\\Extended.seq.context.21bases.pdf"))
  expect_equal(out, TRUE)
  expect_equal(out1, TRUE)
  unlink(paste0(tempdir(), "\\Extended.seq.context.21bases.pdf"))
})

test_that("extended seq.context functions for Strelka GRCh38 vcf", {
  vcf <- ReadAndSplitStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh38.vcf")
  sbs.vcf <- vcf$SBS[[1]]
  mat <- CreateOnePPMFromSBSVCF(sbs.vcf, ref.genome = "GRCh38",
                                seq.context.width = 10)
  out <- PlotPPM(mat, title = "ExtendedSeqContext_21bases")
  out1 <- 
    PlotPPMToPdf(list(mat), titles = "ExtendedSeqContext_21bases",
                 file = paste0(tempdir(), "\\Extended.seq.context.21bases.pdf"))
  expect_equal(out, TRUE)
  expect_equal(out1, TRUE)
  unlink(paste0(tempdir(), "\\Extended.seq.context.21bases.pdf"))
})

test_that("extended seq.context functions for Strelka GRCh37 vcf", {
  vcf <- ReadAndSplitStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh37.vcf")
  sbs.vcf <- vcf$SBS[[1]]
  mat <- CreateOnePPMFromSBSVCF(sbs.vcf, ref.genome = "GRCh37",
                                seq.context.width = 10)
  out <- PlotPPM(mat, title = "ExtendedSeqContext_21bases")
  out1 <- 
    PlotPPMToPdf(list(mat), titles = "ExtendedSeqContext_21bases",
                 file = paste0(tempdir(), "\\Extended.seq.context.21bases.pdf"))
  expect_equal(out, TRUE)
  expect_equal(out1, TRUE)
  unlink(paste0(tempdir(), "\\Extended.seq.context.21bases.pdf"))
})
