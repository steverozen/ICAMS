context("PlotCatalogToPdf.DBS136Catalog")

test_that("PlotCatalogToPdf.DBS136Catalog function is working properly", {
  catalog.counts <- ReadCatalog("testdata/regress.cat.dbs.136.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)
  out <-
    PlotCatalogToPdf(catalog.counts, 
                     file = paste0(tempdir(), "\\PlotCatDBS136.counts.test.pdf"))
  expect_equal(out, TRUE)

  catalog.density <-
    TransformCatalog(catalog.counts, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density")
  out <-
    PlotCatalogToPdf(catalog.density, 
                     file = paste0(tempdir(), "\\PlotCatDBS136.density.test.pdf"))
  expect_equal(out, TRUE)

  unlink(paste0(tempdir(), "\\PlotCatDBS136.counts.test.pdf"))
  unlink(paste0(tempdir(), "\\PlotCatDBS136.density.test.pdf"))
  graphics.off()
})
