context("extended seq.context functions")

test_that("extended seq.context functions are working properly for Mutect GRCh38 vcf", {
  vcf <- ReadAndSplitMutectVCFs("testdata/Mutect.GRCh38.vcf")
  sbs.vcf <- vcf$SBS[[1]]
  mat <- CreateOnePWMFromSBSVCF(sbs.vcf, ref.genome = "GRCh38",
                                seq.context.width = 10)
  out <- PlotPWM(mat, title = "ExtendedSeqContext_21bases")
  out1 <- PlotPWMToPdf(list(mat), titles = "ExtendedSeqContext_21bases",
                       file = "Extended.seq.context.21bases.pdf")
  expect_equal(out, TRUE)
  expect_equal(out1, TRUE)
  unlink("Extended.seq.context.21bases.pdf")
})

test_that("extended seq.context functions are working properly for Mutect GRCh37 vcf", {
  vcf <- ReadAndSplitMutectVCFs("testdata/Mutect.GRCh37.vcf")
  sbs.vcf <- vcf$SBS[[1]]
  mat <- CreateOnePWMFromSBSVCF(sbs.vcf, ref.genome = "GRCh37",
                                seq.context.width = 10)
  out <- PlotPWM(mat, title = "ExtendedSeqContext_21bases")
  out1 <- PlotPWMToPdf(list(mat), titles = "ExtendedSeqContext_21bases",
                       file = "Extended.seq.context.21bases.pdf")
  expect_equal(out, TRUE)
  expect_equal(out1, TRUE)
  unlink("Extended.seq.context.21bases.pdf")
})

test_that("extended seq.context functions are working properly for Strelka GRCh38 vcf", {
  vcf <- ReadAndSplitStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh38.vcf")
  sbs.vcf <- vcf$SBS[[1]]
  mat <- CreateOnePWMFromSBSVCF(sbs.vcf, ref.genome = "GRCh38",
                                seq.context.width = 10)
  out <- PlotPWM(mat, title = "ExtendedSeqContext_21bases")
  out1 <- PlotPWMToPdf(list(mat), titles = "ExtendedSeqContext_21bases",
                       file = "Extended.seq.context.21bases.pdf")
  expect_equal(out, TRUE)
  expect_equal(out1, TRUE)
  unlink("Extended.seq.context.21bases.pdf")
})

test_that("extended seq.context functions are working properly for Strelka GRCh37 vcf", {
  vcf <- ReadAndSplitStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh37.vcf")
  sbs.vcf <- vcf$SBS[[1]]
  mat <- CreateOnePWMFromSBSVCF(sbs.vcf, ref.genome = "GRCh37",
                                seq.context.width = 10)
  out <- PlotPWM(mat, title = "ExtendedSeqContext_21bases")
  out1 <- PlotPWMToPdf(list(mat), titles = "ExtendedSeqContext_21bases",
                       file = "Extended.seq.context.21bases.pdf")
  expect_equal(out, TRUE)
  expect_equal(out1, TRUE)
  unlink("Extended.seq.context.21bases.pdf")
})