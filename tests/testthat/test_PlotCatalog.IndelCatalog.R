context("PlotCatalog.IndelCatalog")

test_that("PlotCatalog.IndelCatalog function", {
  opar <- par(no.readonly = TRUE)
  par(mar = c(2, 2, 2, 1))
  catalog <- ReadCatalog("testdata/BTSG_WGS_PCAWG.indels.csv",
                       ref.genome = "GRCh37",
                       region = "genome", catalog.type = "counts")
  cat <- catalog[, 1, drop = FALSE]
  out <- PlotCatalog(cat)
  expect_equal(out, TRUE)

  cat.counts.signature <- apply(cat, MARGIN = 2, function(x) x / sum(x))
  cat.counts.signature <-
    as.catalog(cat.counts.signature, ref.genome = "GRCh37",
               region = "genome", catalog.type = "counts.signature")
  out <- PlotCatalog(cat.counts.signature)
  expect_equal(out, TRUE)
  on.exit(par(opar))
})
