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
  expect_equal(out$plot.success, TRUE)

  cat.counts.signature <- apply(cat, MARGIN = 2, function(x) x / sum(x))
  cat.counts.signature <-
    as.catalog(cat.counts.signature, ref.genome = "GRCh37",
               region = "genome", catalog.type = "counts.signature")
  out <- PlotCatalog(cat.counts.signature)
  expect_equal(out$plot.success, TRUE)
  graphics.off()
})
