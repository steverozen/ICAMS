context("PlotCatalog.DBS78Catalog")

test_that("PlotCatalog.DBS78Catalog function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  opar <- par(mar = c(2, 2, 2, 1))
  on.exit(par(opar))
  catalog.counts <- ReadCatalog("testdata/regress.cat.dbs.78.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  cat.counts <- catalog.counts[, 1, drop = FALSE]
  out <- PlotCatalog(cat.counts)
  out1 <- PlotCatalog(cat.counts, grid = FALSE, upper = FALSE, xlabels = FALSE,
                     ylabels = FALSE)
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)

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
  out <- PlotCatalog(cat.counts.signature)
  colnames(cat.counts.signature) <- ""
  out1 <- PlotCatalog(cat.counts.signature, grid = FALSE, upper = FALSE, xlabels = FALSE,
                      ylabels = FALSE)
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)

  cat.density.signature <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")
  out <- PlotCatalog(cat.density.signature)
  expect_equal(out$plot.success, TRUE)
  graphics.off()
})
