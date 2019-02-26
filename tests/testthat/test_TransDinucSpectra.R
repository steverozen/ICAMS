context("TransDinucSpectra")

test_that("TransDinucSpectra function is working properly", {
  human.catDNS78 <- ReadCatDNS78("testdata/regress.cat.dns.78.csv")
  mouse.catDNS78 <- ReadCatDNS78("testdata/mouse.cat.dns.78.csv")
  inferred.mouse.catDNS78 <-
    TransDinucSpectra(human.catDNS78,
                      source.abundance = "GRCh37.genome",
                      target.abundance = "GRCm38.genome")
  expect_equal(inferred.mouse.catDNS78, mouse.catDNS78)
})
