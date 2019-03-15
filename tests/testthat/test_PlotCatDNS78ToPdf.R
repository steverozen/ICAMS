context("PlotCatDNS78ToPdf")

test_that("PlotCatDNS78ToPdf function is working properly", {
  catalog <- ReadCatDNS78("testdata/regress.cat.dns.78.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  cat1 <- catalog[, 1, drop = FALSE]
  cat1.sig <-
    TransformCatalog(cat1, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome", which.n = 2,
                     source.type = "counts", target.type = "signature")
  cat1.density <-
    TransformCatalog(cat1, source.abundance = "GRCh37.genome",
                     source.type = "counts", target.type = "density", which.n = 2)

  cat2 <- catalog[, 2, drop = FALSE]
  cat2.sig <-
    TransformCatalog(cat2, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome", which.n = 2,
                     source.type = "counts", target.type = "signature")
  cat2.density <-
    TransformCatalog(cat2, source.abundance = "GRCh37.genome",
                     source.type = "counts", target.type = "density", which.n = 2)

  cat3 <- catalog[, 3, drop = FALSE]
  cat3.sig <-
    TransformCatalog(cat3, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome", which.n = 2,
                     source.type = "counts", target.type = "signature")
  cat3.density <-
    TransformCatalog(cat3, source.abundance = "GRCh37.genome",
                     source.type = "counts", target.type = "density", which.n = 2)

  cat4 <- catalog[, 4, drop = FALSE]
  cat4.sig <-
    TransformCatalog(cat4, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome", which.n = 2,
                     source.type = "counts", target.type = "signature")
  cat4.density <-
    TransformCatalog(cat4, source.abundance = "GRCh37.genome",
                     source.type = "counts", target.type = "density", which.n = 2)

  cat <- cbind(cat1, cat1.sig, cat1.density,
               cat2, cat2.sig, cat2.density,
               cat3, cat3.sig, cat3.density,
               cat4, cat4.sig, cat4.density)

  type <- c("counts", "signature", "density")
  out <-
    PlotCatDNS78ToPdf(cat, "PlotCatDNS78.test.pdf",
                      type = rep(type, 4))
  expect_equal(out, TRUE)
  unlink("PlotCatDNS78.test.pdf")
  graphics.off()
})
