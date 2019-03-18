context("PlotSNSClassStrandBiasToPdf")

test_that("PlotSNSClassStrandBiasToPdf function is working properly", {
  catalog <- ReadCatSNS192("testdata/regress.cat.sns.192.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  cat1 <- catalog[, 1, drop = FALSE]
  cat1.sig <-
    TransformCatalog(cat1, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome", which.n = 3,
                     source.type = "counts", target.type = "signature")
  cat2 <- catalog[, 2, drop = FALSE]
  cat2.sig <-
    TransformCatalog(cat2, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome", which.n = 3,
                     source.type = "counts", target.type = "signature")
  cat3 <- catalog[, 3, drop = FALSE]
  cat3.sig <-
    TransformCatalog(cat3, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome", which.n = 3,
                     source.type = "counts", target.type = "signature")
  cat4 <- catalog[, 4, drop = FALSE]
  cat4.sig <-
    TransformCatalog(cat4, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome", which.n = 3,
                     source.type = "counts", target.type = "signature")

  cat <- cbind(cat1, cat1.sig,
               cat2, cat2.sig,
               cat3, cat3.sig,
               cat4, cat4.sig)

  type <- c("counts", "signature")
  out <-
    PlotSNSClassStrandBiasToPdf(cat,
                                filename = "PlotSNSClassStrandBias.test.pdf",
                                type = rep(type, 4))
  expect_equal(out, TRUE)
  unlink("PlotSNSClassStrandBias.test.pdf")
  graphics.off()
})
