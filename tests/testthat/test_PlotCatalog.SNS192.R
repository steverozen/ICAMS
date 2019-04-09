context("PlotCatalog.SNS192")

test_that("PlotCatalog.SNS192 function is working properly", {
  par(mar = c(6, 4, 4, 1))
  catalog.counts <-
    ReadCatalog("testdata/regress.cat.sns.192.csv", ref.genome = "GRCh37",
                         region = "transcript", catalog.type = "counts")
  cat.counts <- catalog.counts[, 1, drop = FALSE]
  cat.counts <- PreserveCatalogAttribute(catalog.counts, cat.counts)
  out <- PlotCatalog(cat.counts)
  expect_equal(out, TRUE)

  cat.density <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "transcript",
                     target.catalog.type = "density")
  out <- PlotCatalog(cat.density)
  expect_equal(out, TRUE)

  cat.counts.signature <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "transcript",
                     target.catalog.type = "counts.signature")
  out <- PlotCatalog(cat.counts.signature)
  expect_equal(out, TRUE)

  cat.density.signature <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "transcript",
                     target.catalog.type = "density.signature")
  out <- PlotCatalog(cat.density.signature)
  expect_equal(out, TRUE)

})
