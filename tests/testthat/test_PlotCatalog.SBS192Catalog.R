context("PlotCatalog.SBS192Catalog")

test_that("PlotCatalog.SBS192Catalog function is working properly", {
  par(mar = c(2, 2, 2, 1))
  catalog.counts <-
    ReadCatalog("testdata/regress.cat.sbs.192.csv", 
                ref.genome = "GRCh37",
                region = "transcript", catalog.type = "counts")
  cat.counts <- catalog.counts[, 1, drop = FALSE]
  out <- PlotCatalog(cat.counts)
  out1 <- PlotCatalog(cat.counts, plot.SBS12 = TRUE)
  expect_equal(out, TRUE)
  expect_equal(out1, TRUE)

  cat.density <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "transcript",
                     target.catalog.type = "density")
  out <- PlotCatalog(cat.density)
  out1 <- PlotCatalog(cat.density, plot.SBS12 = TRUE)
  expect_equal(out, TRUE)
  expect_equal(out1, TRUE)

  cat.counts.signature <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "transcript",
                     target.catalog.type = "counts.signature")
  out <- PlotCatalog(cat.counts.signature)
  out1 <- PlotCatalog(cat.counts.signature, plot.SBS12 = TRUE)
  expect_equal(out, TRUE)
  expect_equal(out1, TRUE)

  cat.density.signature <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "transcript",
                     target.catalog.type = "density.signature")
  out <- PlotCatalog(cat.density.signature)
  out1 <- PlotCatalog(cat.density.signature, plot.SBS12 = TRUE)
  expect_equal(out, TRUE)
  expect_equal(out1, TRUE)
  })
