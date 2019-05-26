context("PlotCatalog.SNS192CatalogNoContext")

test_that("PlotCatalog.SNS192CatalogNoContext function is working properly", {
  par(mar = c(5, 5, 5, 1))
  catalog.counts <- ReadCatalog("testdata/regress.cat.sns.192.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  cat.counts <- catalog.counts[, 1, drop = FALSE]
  out <- PlotCatalog(cat.counts, no.context = TRUE)
  expect_equal(out, TRUE)

  cat.density <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density")
  out <- PlotCatalog(cat.density, no.context = TRUE)
  expect_equal(out, TRUE)

  cat.counts.signature <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")
  out <- PlotCatalog(cat.counts.signature, no.context = TRUE)
  expect_equal(out, TRUE)

  cat.density.signature <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")
  out <- PlotCatalog(cat.density.signature, no.context = TRUE)
  expect_equal(out, TRUE)
})
