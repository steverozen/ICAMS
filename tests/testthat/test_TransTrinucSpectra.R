context("TransTrinucSpectra")

test_that("TransTrinucSpectra function is working properly", {
  human.catSNS96 <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")
  mouse.catSNS96 <- ReadCatSNS96("testdata/mouse.cat.sns.96.csv")
  inferred.mouse.catSNS96 <- TransTrinucSpectra(human.catSNS96,
                                                source.abundance = "GRCh37.genome",
                                                target.abundance = "GRCm38.genome")
  expect_equal(inferred.mouse.catSNS96, mouse.catSNS96)
})
