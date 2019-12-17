context("VCFFilestoZipfile functions")

test_that("MutectVCFFilesToZipFile function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  file <- "testdata/Mutect-GRCh37/"
  out <- MutectVCFFilesToZipFile(file, ref.genome = "hg19",
                                 trans.ranges = trans.ranges.GRCh37,
                                 region = "genome",
                                 zipfile.name = "test")
  expect_type(out, "list")
  name <- grep(".zip", list.files("testdata/Mutect-GRCh37/"), value = TRUE)
  expect_equal(name, "test.zip")
  unlink("testdata/Mutect-GRCh37/test.zip")
  graphics.off()
  unlink("testdata/Mutect-GRCh37/Rplots.pdf")
})

test_that("StrelkaSBSVCFFilesToZipFile function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  file <- "testdata/Strelka-SBS-GRCh37/"
  out <- StrelkaSBSVCFFilesToZipFile(file, ref.genome = "hg19",
                                     trans.ranges = trans.ranges.GRCh37,
                                     region = "genome",
                                     zipfile.name = "test")
  expect_type(out, "list")
  name <- grep(".zip", list.files("testdata/Strelka-SBS-GRCh37/"), value = TRUE)
  expect_equal(name, "test.zip")
  unlink("testdata/Strelka-SBS-GRCh37/test.zip")
  graphics.off()
  unlink("testdata/Strelka-SBS-GRCh37/Rplots.pdf")
})

test_that("StrelkaIDVCFFilesToZipFile function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  file <- "testdata/Strelka-ID-GRCh37/"
  out <- StrelkaIDVCFFilesToZipFile(file, ref.genome = "hg19",
                                    region = "genome",
                                    zipfile.name = "test")
  expect_type(out, "list")
  name <- grep(".zip", list.files("testdata/Strelka-ID-GRCh37/"), value = TRUE)
  expect_equal(name, "test.zip")
  unlink("testdata/Strelka-ID-GRCh37/test.zip")
  graphics.off()
  unlink("testdata/Strelka-ID-GRCh37/Rplots.pdf")
})