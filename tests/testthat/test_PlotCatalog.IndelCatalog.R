context("PlotCatalog.IndelCatalog")

test_that("PlotCatalog.IndelCatalog function", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  opar <- par(mar = c(2, 2, 2, 1))
  on.exit(par(opar))
  catalog <- ReadCatalog("testdata/BTSG_WGS_PCAWG.indels.csv",
                       ref.genome = "GRCh37",
                       region = "genome", catalog.type = "counts")
  cat <- catalog[, 1, drop = FALSE]
  out <- PlotCatalog(cat)
  out1 <- PlotCatalog(cat, ylim = c(0, 1000))
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)

  cat.counts.signature <- apply(cat, MARGIN = 2, function(x) x / sum(x))
  cat.counts.signature <-
    as.catalog(cat.counts.signature, ref.genome = "GRCh37",
               region = "genome", catalog.type = "counts.signature")
  out <- PlotCatalog(cat.counts.signature)
  out1 <- PlotCatalog(cat.counts.signature, ylim = c(0, 0.5))
  expect_equal(out$plot.success, TRUE)
  expect_equal(out1$plot.success, TRUE)
  graphics.off()
})
