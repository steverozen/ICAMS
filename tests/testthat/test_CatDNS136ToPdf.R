context("CatDNS136ToPdf")

test_that("CatDNS136ToPdf function is working properly", {
  catalog <- ReadCatDNS136("testdata/regress.cat.dns.136.csv")
  colnames(catalog) <- paste0("HepG2_", 1 : 4)
  type <- rep(c("density", "counts"), 2)
  out <-
    CatDNS136ToPdf(catalog, "PlotCatDNS136.test.pdf",
                   type = type,
                   abundance = abundance.4bp.genome.GRCh37)
  expect_equal(out, TRUE)
})
