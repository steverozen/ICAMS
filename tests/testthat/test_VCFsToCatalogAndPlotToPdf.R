context("VCFsToCatalogsAndPlotToPdf functions")

test_that("CatalogAndPlotToPdf function for Mutect VCFs", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  files <- list.files(path = "testdata/Mutect-GRCh37", full.names = TRUE)
  catalogs1 <- VCFsToCatalogsAndPlotToPdf(files, ref.genome = "hg19", 
                                         output.dir = tempdir(),
                                         variant.caller = "mutect",
                                         region = "genome",
                                         base.filename = "mutect")
  names <- grep("\\.pdf$", list.files(tempdir()), value = TRUE)
  expect_equal(length(names), 9)
  unlink(file.path(tempdir(), names))
})

test_that("CatalogAndPlotToPdf function for Strelka SBS VCFs", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  files <- list.files(path = "testdata/Strelka-SBS-GRCh37/", full.names = TRUE)
  catalogs1 <- VCFsToCatalogsAndPlotToPdf(files, ref.genome = "hg19", 
                                         output.dir = tempdir(),
                                         variant.caller = "strelka",
                                         region = "genome",
                                         base.filename = "strelka.sbs")
  names <- grep("\\.pdf$", list.files(tempdir()), value = TRUE)
  expect_equal(length(names), 7)
  unlink(file.path(tempdir(), names))
})

test_that("CatalogAndPlotToPdf function for Strelka ID VCFs", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  files <- list.files(path = "testdata/Strelka-ID-GRCh37/", full.names = TRUE)
  catalogs1 <- VCFsToCatalogsAndPlotToPdf(files, ref.genome = "hg19", 
                                         output.dir = tempdir(),
                                         variant.caller = "strelka",
                                         region = "genome",
                                         base.filename = "strelka.id")
  names <- grep("\\.pdf$", list.files(tempdir()), value = TRUE)
  expect_equal(length(names), 2)
  unlink(file.path(tempdir(), names))
})