context("PlotCatalog.DBS144Catalog")

test_that("PlotCatalog.DBS144Catalog function", {
  opar <- par(no.readonly = TRUE)
  par(mar = c(1, 1, 1, 1))
  catalog.counts <-
    ReadCatalog("testdata/regress.cat.dbs.144.csv",
                ref.genome = "GRCh37",
                region = "transcript",
                catalog.type = "counts")
  cat.counts <- catalog.counts[, 1, drop = FALSE]
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
  on.exit(par(opar))
})
