context("[ methods for catalogs")

test_that("[ method for SNS96Catalog is working properly", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.sns.96.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out1 <- PlotCatalog(catalog.counts[, 1, drop = FALSE])
  out2 <- PlotCatalogToPdf(catalog.counts[, 1, drop = FALSE],
                           file = "test.pdf")
  expect_equal(out1, TRUE)
  expect_equal(out2, TRUE)
  expect_error(PlotCatalog(catalog.counts[, 1]))
  unlink("test.pdf")
})

test_that("[ method for SNS192Catalog is working properly", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.sns.192.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out1 <- PlotCatalog(catalog.counts[, 1, drop = FALSE])
  out2 <- PlotCatalogToPdf(catalog.counts[, 1, drop = FALSE],
                           file = "test.pdf")
  expect_equal(out1, TRUE)
  expect_equal(out2, TRUE)
  expect_error(PlotCatalog(catalog.counts[, 1]))
  unlink("test.pdf")
})

test_that("[ method for SNS1536Catalog is working properly", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.sns.1536.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out1 <- PlotCatalog(catalog.counts[, 1, drop = FALSE])
  out2 <- PlotCatalogToPdf(catalog.counts[, 1, drop = FALSE],
                           file = "test.pdf")
  expect_equal(out1, TRUE)
  expect_equal(out2, TRUE)
  expect_error(PlotCatalog(catalog.counts[, 1]))
  unlink("test.pdf")
})

test_that("[ method for DNS78Catalog is working properly", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.dns.78.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out1 <- PlotCatalog(catalog.counts[, 1, drop = FALSE])
  out2 <- PlotCatalogToPdf(catalog.counts[, 1, drop = FALSE],
                           file = "test.pdf")
  expect_equal(out1, TRUE)
  expect_equal(out2, TRUE)
  expect_error(PlotCatalog(catalog.counts[, 1]))
  unlink("test.pdf")
})

test_that("[ method for DNS144Catalog is working properly", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.dns.144.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out1 <- PlotCatalog(catalog.counts[, 1, drop = FALSE])
  out2 <- PlotCatalogToPdf(catalog.counts[, 1, drop = FALSE],
                           file = "test.pdf")
  expect_equal(out1, TRUE)
  expect_equal(out2, TRUE)
  expect_error(PlotCatalog(catalog.counts[, 1]))
  unlink("test.pdf")
})

test_that("[ method for DNS136Catalog is working properly", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.dns.136.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out1 <- PlotCatalog(catalog.counts[, 1, drop = FALSE])
  out2 <- PlotCatalogToPdf(catalog.counts[, 1, drop = FALSE],
                           file = "test.pdf")
  expect_equal(out1, TRUE)
  expect_equal(out2, TRUE)
  expect_error(PlotCatalog(catalog.counts[, 1]))
  unlink("test.pdf")
})

test_that("[ method for IndelCatalog is working properly", {
  catalog.counts <- ReadCatalog("testdata/BTSG_WGS_PCAWG.indels.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  out1 <- PlotCatalog(catalog.counts[, 1, drop = FALSE])
  out2 <- PlotCatalogToPdf(catalog.counts[, 1, drop = FALSE],
                           file = "test.pdf")
  expect_equal(out1, TRUE)
  expect_equal(out2, TRUE)
  expect_error(PlotCatalog(catalog.counts[, 1]))
  unlink("test.pdf")
})
