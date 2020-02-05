context("PlotCatalog.SBS96Catalog")

test_that("PlotCatalog.SBS96Catalog for one column catalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  opar <- par(mar = c(5.5, 5, 5, 1))
  on.exit(par(opar))
  catalog.counts <-
    ReadCatalog("testdata/regress.cat.sbs.96.csv", ref.genome = "GRCh37",
                region = "genome", catalog.type = "counts")
  cat.counts <- catalog.counts[, 1, drop = FALSE]
  out <- PlotCatalog(cat.counts)
  out1 <- PlotCatalog(cat.counts, ylim = c(0, 1500))
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)

  cat.density <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density")
  out <- PlotCatalog(cat.density)
  out1 <- PlotCatalog(cat.density, ylim = c(0, 20))
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)

  cat.counts.signature <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")
  out <- PlotCatalog(cat.counts.signature)
  out1 <- PlotCatalog(cat.counts.signature, ylim = c(0, 0.2))
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)

  cat.density.signature <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")
  out <- PlotCatalog(cat.density.signature)
  out1 <- PlotCatalog(cat.density.signature, ylim = c(0, 0.2))
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)
  graphics.off()
})

test_that("PlotCatalog.SBS96Catalog for two column catalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  opar <- par(mar = c(5.5, 5, 5, 1))
  on.exit(par(opar))
  catalog.counts <-
    ReadCatalog("testdata/regress.cat.sbs.96.csv", ref.genome = "GRCh37",
                region = "genome", catalog.type = "counts")
  cat.counts <- catalog.counts[, 1:2]
  out <- PlotCatalog(cat.counts)
  out1 <- PlotCatalog(cat.counts, ylim = c(0, 2500))
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)
  
  cat.density <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density")
  out <- PlotCatalog(cat.density)
  out1 <- PlotCatalog(cat.density, ylim = c(0, 35))
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)
  
  cat.counts.signature <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")
  out <- PlotCatalog(cat.counts.signature)
  out1 <- PlotCatalog(cat.counts.signature, ylim = c(0, 0.3))
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)
  
  cat.density.signature <-
    TransformCatalog(cat.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")
  out <- PlotCatalog(cat.density.signature)
  out1 <- PlotCatalog(cat.density.signature, ylim = c(0, 0.3))
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)
  graphics.off()
})