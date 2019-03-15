context("PlotCatDNS136ToPdf")

test_that("PlotCatDNS136ToPdf function is working properly", {
  catalog <- ReadCatDNS136("testdata/regress.cat.dns.136.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  cat1 <- catalog[, 1, drop = FALSE]
  cat1.density <-
    TransformCatalog(cat1, source.abundance = "GRCh37.genome",
                     which.n = 4,
                     source.type = "counts", target.type = "density")
  cat2 <- catalog[, 2, drop = FALSE]
  cat2.density <-
    TransformCatalog(cat2, source.abundance = "GRCh37.genome",
                     which.n = 4,
                     source.type = "counts", target.type = "density")
  cat3 <- catalog[, 3, drop = FALSE]
  cat3.density <-
    TransformCatalog(cat3, source.abundance = "GRCh37.genome",
                     which.n = 4,
                     source.type = "counts", target.type = "density")
  cat4 <- catalog[, 4, drop = FALSE]
  cat4.density <-
    TransformCatalog(cat4, source.abundance = "GRCh37.genome",
                     which.n = 4,
                     source.type = "counts", target.type = "density")

  cat <- cbind(cat1, cat1.density,
               cat2, cat2.density,
               cat3, cat3.density,
               cat4, cat4.density)

  type <- c("counts", "density")
  out <-
    PlotCatDNS136ToPdf(cat, name = "PlotCatDNS136.test.pdf",
                       type = rep(type, 4))
  expect_equal(out, TRUE)
  unlink("PlotCatDNS136.test.pdf")
  graphics.off()
})
