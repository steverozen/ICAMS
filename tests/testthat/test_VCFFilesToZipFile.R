context("VCFFilestoZipfile functions")

test_that("MutectVCFFilesToZipFile function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  dir <- "testdata/Mutect-GRCh37"
  out <- MutectVCFFilesToZipFile(dir, 
                                 file = tempdir(), 
                                 ref.genome = "hg19",
                                 trans.ranges = trans.ranges.GRCh37,
                                 region = "genome",
                                 zipfile.name = "test")
  expect_type(out, "list")
  name <- grep(".zip", list.files(tempdir()), value = TRUE)
  expect_equal(name, "test.zip")
  unlink(paste0(tempdir(), "/test.zip"))
  graphics.off()
  unlink("testdata/Mutect-GRCh37/Rplots.pdf")
})

test_that("StrelkaSBSVCFFilesToZipFile function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  dir <- "testdata/Strelka-SBS-GRCh37"
  out <- StrelkaSBSVCFFilesToZipFile(dir,
                                     zipfile = paste0(tempdir(), "/test.zip"), 
                                     ref.genome = "hg19",
                                     trans.ranges = trans.ranges.GRCh37,
                                     region = "genome")
  expect_type(out, "list")
  name <- grep("\\.zip$", list.files(tempdir()), value = TRUE)
  expect_equal(name, "test.zip")
  unlink(paste0(tempdir(), "/test.zip"))
  graphics.off()
  unlink("testdata/Strelka-SBS-GRCh37/Rplots.pdf")
})

test_that("StrelkaIDVCFFilesToZipFile function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  dir <- "testdata/Strelka-ID-GRCh37"
  out <- StrelkaIDVCFFilesToZipFile(dir,
                                    file = tempdir(),
                                    ref.genome = "hg19",
                                    region = "genome",
                                    zipfile.name = "test")
  expect_type(out, "list")
  name <- grep(".zip", list.files(tempdir()), value = TRUE)
  expect_equal(name, "test.zip")
  unlink(paste0(tempdir(), "/test.zip"))
  graphics.off()
  unlink("testdata/Strelka-ID-GRCh37/Rplots.pdf")
})