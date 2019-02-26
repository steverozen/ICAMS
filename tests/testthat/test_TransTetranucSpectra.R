context("TransTetranucSpectra")

test_that("TransTetranucSpectra function is working properly", {
  human.catDNS136 <- ReadCatDNS136("testdata/regress.cat.dns.136.csv")
  mouse.catDNS136 <- ReadCatDNS136("testdata/mouse.cat.dns.136.csv")
  inferred.mouse.catDNS136 <-
    TransTetranucSpectra(human.catDNS136,
                         source.abundance = "GRCh37.genome",
                         target.abundance = "GRCm38.genome")
  expect_equal(inferred.mouse.catDNS136, mouse.catDNS136)
})
