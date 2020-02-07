context("PlotCatalog.DBS144Catalog")

test_that("PlotCatalog.DBS144Catalog function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  opar <- par(mar = c(2, 4, 2, 2))
  on.exit(par(opar))
  catalog.counts <-
    ReadCatalog("testdata/regress.cat.dbs.144.csv",
                ref.genome = "GRCh37",
                region = "transcript",
                catalog.type = "counts")
  cat.counts <- catalog.counts[, 1, drop = FALSE]
  out <- PlotCatalog(cat.counts)
  out1 <- PlotCatalog(cat.counts, cex = 0.8)
  graphics.off()
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)
  
  cat.density <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "transcript",
                     target.catalog.type = "density")
  out <- PlotCatalog(cat.density)
  out1 <- PlotCatalog(cat.density, cex = 0.8)
  graphics.off()
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)

  cat.counts.signature <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "transcript",
                     target.catalog.type = "counts.signature")
  out <- PlotCatalog(cat.counts.signature)
  out1 <- PlotCatalog(cat.counts.signature, cex = 0.8)
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)

  cat.density.signature <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "transcript",
                     target.catalog.type = "density.signature")
  out <- PlotCatalog(cat.density.signature)
  out1 <- PlotCatalog(cat.density.signature, cex = 0.8)
  graphics.off()
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)
})
