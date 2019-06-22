context("[ methods for catalogs")

test_that("[ method for SBS96Catalog", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out1 <- PlotCatalog(catalog.counts[, 1, drop = FALSE])
  out2 <- PlotCatalogToPdf(catalog.counts[, 1, drop = FALSE],
                           file = file.path(tempdir(), "test.pdf"))
  expect_equal(out1, TRUE)
  expect_equal(out2, TRUE)
  expect_error(PlotCatalog(catalog.counts[, 1]))
  unlink(file.path(tempdir(), "test.pdf"))
})

test_that("[ method for SBS192Catalog", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.sbs.192.csv",
                                ref.genome = "GRCh37",
                                region = "transcript",
                                catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out1 <- PlotCatalog(catalog.counts[, 1, drop = FALSE])
  out2 <- PlotCatalogToPdf(catalog.counts[, 1, drop = FALSE],
                           file = file.path(tempdir(), "test.pdf"))
  expect_equal(out1, TRUE)
  expect_equal(out2, TRUE)
  expect_error(PlotCatalog(catalog.counts[, 1]))
  unlink(file.path(tempdir(), "test.pdf"))
})

test_that("[ method for SBS1536Catalog", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.sbs.1536.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out1 <- PlotCatalog(catalog.counts[, 1, drop = FALSE])
  out2 <- PlotCatalogToPdf(catalog.counts[, 1, drop = FALSE],
                           file = file.path(tempdir(), "test.pdf"))
  expect_equal(out1, TRUE)
  expect_equal(out2, TRUE)
  expect_error(PlotCatalog(catalog.counts[, 1]))
  unlink(file.path(tempdir(), "test.pdf"))
})

test_that("[ method for DBS78Catalog", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.dbs.78.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out1 <- PlotCatalog(catalog.counts[, 1, drop = FALSE])
  out2 <- PlotCatalogToPdf(catalog.counts[, 1, drop = FALSE],
                           file = file.path(tempdir(), "test.pdf"))
  expect_equal(out1, TRUE)
  expect_equal(out2, TRUE)
  expect_error(PlotCatalog(catalog.counts[, 1]))
  unlink(file.path(tempdir(), "test.pdf"))
})

test_that("[ method for DBS144Catalog", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.dbs.144.csv",
                                ref.genome = "GRCh37",
                                region = "transcript",
                                catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out1 <- PlotCatalog(catalog.counts[, 1, drop = FALSE])
  out2 <- PlotCatalogToPdf(catalog.counts[, 1, drop = FALSE],
                           file = file.path(tempdir(), "test.pdf"))
  expect_equal(out1, TRUE)
  expect_equal(out2, TRUE)
  expect_error(PlotCatalog(catalog.counts[, 1]))
  unlink(file.path(tempdir(), "test.pdf"))
})

test_that("[ method for DBS136Catalog", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.dbs.136.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out1 <- PlotCatalog(catalog.counts[, 1, drop = FALSE])
  out2 <- PlotCatalogToPdf(catalog.counts[, 1, drop = FALSE],
                           file = file.path(tempdir(), "test.pdf"))
  expect_equal(out1, TRUE)
  expect_equal(out2, TRUE)
  expect_error(PlotCatalog(catalog.counts[, 1]))
  unlink(file.path(tempdir(), "test.pdf"))
})

test_that("[ method for IndelCatalog", {
  catalog.counts <- ReadCatalog("testdata/BTSG_WGS_PCAWG.indels.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  out1 <- PlotCatalog(catalog.counts[, 1, drop = FALSE])
  out2 <- PlotCatalogToPdf(catalog.counts[, 1, drop = FALSE],
                           file = file.path(tempdir(), "test.pdf"))
  expect_equal(out1, TRUE)
  expect_equal(out2, TRUE)
  expect_error(PlotCatalog(catalog.counts[, 1]))
  unlink(file.path(tempdir(), "test.pdf"))
})

