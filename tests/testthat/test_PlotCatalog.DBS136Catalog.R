context("PlotCatalog.DBS136Catalog")

test_that("PlotCatalog.DBS136Catalog function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  opar <- par(oma = c(2, 2, 2, 0))
  on.exit(par(opar))
  # cat("just after on.exit, oma =", par("oma"), "\n")
  catalog.counts <- ReadCatalog("testdata/regress.cat.dbs.136.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  cat.counts <- catalog.counts[, 1, drop = FALSE]
  out <- PlotCatalog(cat.counts)
  graphics.off()
  expect_equal(out$plot.success, TRUE)

  cat.density <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density")
  out <- PlotCatalog(cat.density)
  expect_equal(out$plot.success, TRUE)

  cat.counts.signature <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")
  expect_error(PlotCatalog(cat.counts.signature))

  cat.density.signature <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")
  expect_error(PlotCatalog(cat.density.signature))
  graphics.off()
})

# cat("after testthat, oma = ", par("oma"), "\n")
