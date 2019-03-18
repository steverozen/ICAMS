context("PlotCatSNS96ToPdf")

test_that("PlotCatSNS96ToPdf function is working properly", {
  catalog <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  cat1 <- catalog[, 1, drop = FALSE]
  cat1.sig <-
    TransformCatalog(cat1, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome", which.n = 3,
                     source.type = "counts", target.type = "signature")
  cat1.density <-
    TransformCatalog(cat1, source.abundance = "GRCh37.genome",
                     source.type = "counts", target.type = "density", which.n = 3)

  cat2 <- catalog[, 2, drop = FALSE]
  cat2.sig <-
    TransformCatalog(cat2, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome", which.n = 3,
                     source.type = "counts", target.type = "signature")
  cat2.density <-
    TransformCatalog(cat2, source.abundance = "GRCh37.genome",
                     source.type = "counts", target.type = "density", which.n = 3)

  cat3 <- catalog[, 3, drop = FALSE]
  cat3.sig <-
    TransformCatalog(cat3, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome", which.n = 3,
                     source.type = "counts", target.type = "signature")
  cat3.density <-
    TransformCatalog(cat3, source.abundance = "GRCh37.genome",
                     source.type = "counts", target.type = "density", which.n = 3)

  cat4 <- catalog[, 4, drop = FALSE]
  cat4.sig <-
    TransformCatalog(cat4, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome", which.n = 3,
                     source.type = "counts", target.type = "signature")
  cat4.density <-
    TransformCatalog(cat4, source.abundance = "GRCh37.genome",
                     source.type = "counts", target.type = "density", which.n = 3)

  cat <- cbind(cat1, cat1.sig, cat1.density,
               cat2, cat2.sig, cat2.density,
               cat3, cat3.sig, cat3.density,
               cat4, cat4.sig, cat4.density)

  type <- c("counts", "signature", "density")
  out <-
    PlotCatSNS96ToPdf(cat, filename = "PlotCatSNS96.test.pdf",
                      type = rep(type, 4), grid = TRUE,
                      upper = TRUE,
                      xlabels = TRUE)
  expect_equal(out, TRUE)
  unlink("PlotCatSNS96.test.pdf")
  graphics.off()
})
