context("PlotCatalog.SNSClassStrandBias")

test_that("PlotCatalog.SNSClassStrandBias function is working properly", {
  par(mar = c(5, 5, 5, 1))
  catalog.counts <- ReadCatalog("testdata/regress.cat.sns.192.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  cat.counts <- catalog.counts[, 1, drop = FALSE]
  cat.counts <- PreserveCatalogAttribute(catalog.counts, cat.counts)
  out <- PlotCatalog(cat.counts, strandbias = TRUE)
  expect_equal(out, TRUE)

  cat.density <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density")
  out <- PlotCatalog(cat.density, strandbias = TRUE)
  expect_equal(out, TRUE)

  cat.counts.signature <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")
  out <- PlotCatalog(cat.counts.signature, strandbias = TRUE)
  expect_equal(out, TRUE)

  cat.density.signature <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")
  out <- PlotCatalog(cat.density.signature, strandbias = TRUE)
  expect_equal(out, TRUE)
})
