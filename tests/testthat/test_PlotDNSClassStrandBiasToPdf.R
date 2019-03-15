context("PlotDNSClassStrandBiasToPdf")

test_that("PlotDNSClassStrandBiasToPdf function is working properly", {
  catalog <- ReadCatDNS144("testdata/regress.cat.dns.144.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  cat1 <- catalog[, 1, drop = FALSE]
  cat1.sig <-
    TransformCatalog(cat1, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome", which.n = 2,
                     source.type = "counts", target.type = "signature")
  cat2 <- catalog[, 2, drop = FALSE]
  cat2.sig <-
    TransformCatalog(cat2, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome", which.n = 2,
                     source.type = "counts", target.type = "signature")
  cat3 <- catalog[, 3, drop = FALSE]
  cat3.sig <-
    TransformCatalog(cat3, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome", which.n = 2,
                     source.type = "counts", target.type = "signature")
  cat4 <- catalog[, 4, drop = FALSE]
  cat4.sig <-
    TransformCatalog(cat4, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome", which.n = 2,
                     source.type = "counts", target.type = "signature")

  cat <- cbind(cat1, cat1.sig,
               cat2, cat2.sig,
               cat3, cat3.sig,
               cat4, cat4.sig)

  type <- c("counts", "signature")
  out <- PlotDNSClassStrandBiasToPdf(cat,
                                     name = "PlotDNSClassStrandBias.test.pdf",
                                     type = rep(type, 4))
  expect_equal(out, TRUE)
  unlink("PlotDNSClassStrandBias.test.pdf")
  graphics.off()
})
