context("cbind methods for catalogs")

test_that("cbind method for SBS96Catalog", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)

  cat1 <- catalog.counts[, 1, drop = FALSE]
  cat2 <- catalog.counts[, 2, drop = FALSE]
  cat3 <- cbind(cat1, cat2)
  out <- PlotCatalogToPdf(cat3, file = file.path(tempdir(), "test.pdf"))
  expect_equal(out, TRUE)
  unlink(file.path(tempdir(), "test.pdf"))
})

test_that("cbind method for SBS192Catalog", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.sbs.192.csv",
                                ref.genome = "GRCh37",
                                region = "transcript",
                                catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)

  cat1 <- catalog.counts[, 1, drop = FALSE]
  cat2 <- catalog.counts[, 2, drop = FALSE]
  cat3 <- cbind(cat1, cat2)
  out <- PlotCatalogToPdf(cat3, file = file.path(tempdir(), "test.pdf"))
  expect_equal(out, TRUE)
  unlink(file.path(tempdir(), "test.pdf"))
})

test_that("cbind method for SBS1536Catalog", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.sbs.1536.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)

  cat1 <- catalog.counts[, 1, drop = FALSE]
  cat2 <- catalog.counts[, 2, drop = FALSE]
  cat3 <- cbind(cat1, cat2)
  out <- PlotCatalogToPdf(cat3, file = file.path(tempdir(), "test.pdf"))
  expect_equal(out, TRUE)
  unlink(file.path(tempdir(), "test.pdf"))
})

test_that("cbind method for DBS78Catalog", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.dbs.78.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)

  cat1 <- catalog.counts[, 1, drop = FALSE]
  cat2 <- catalog.counts[, 2, drop = FALSE]
  cat3 <- cbind(cat1, cat2)
  out <- PlotCatalogToPdf(cat3, file = file.path(tempdir(), "test.pdf"))
  expect_equal(out, TRUE)
  unlink(file.path(tempdir(), "test.pdf"))
})

test_that("cbind method for DBS144Catalog", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.dbs.144.csv",
                                ref.genome = "GRCh37",
                                region = "transcript", 
                                catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)

  cat1 <- catalog.counts[, 1, drop = FALSE]
  cat2 <- catalog.counts[, 2, drop = FALSE]
  cat3 <- cbind(cat1, cat2)
  out <- PlotCatalogToPdf(cat3, file = file.path(tempdir(), "test.pdf"))
  expect_equal(out, TRUE)
  unlink(file.path(tempdir(), "test.pdf"))
})

test_that("cbind method for DBS136Catalog", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.dbs.136.csv",
                                ref.genome = "GRCh37",
                                region = "genome",
                                catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)

  cat1 <- catalog.counts[, 1, drop = FALSE]
  cat2 <- catalog.counts[, 2, drop = FALSE]
  cat3 <- cbind(cat1, cat2)
  out <- PlotCatalogToPdf(cat3, file = file.path(tempdir(), "test.pdf"))
  expect_equal(out, TRUE)
  unlink(file.path(tempdir(), "test.pdf"))
})

test_that("cbind method for IndelCatalog", {
  catalog.counts <- ReadCatalog("testdata/BTSG_WGS_PCAWG.indels.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")

  cat1 <- catalog.counts[, 1, drop = FALSE]
  cat2 <- catalog.counts[, 2, drop = FALSE]
  cat3 <- cbind(cat1, cat2)
  out <- PlotCatalogToPdf(cat3, file = file.path(tempdir(), "test.pdf"))
  expect_equal(out, TRUE)
  unlink(file.path(tempdir(), "test.pdf"))
})

