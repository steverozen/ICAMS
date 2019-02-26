context("TransPentanucSpectra")

test_that("TransPentanucSpectra function is working properly", {
  human.catSNS1536 <- ReadCatSNS1536("testdata/regress.cat.sns.1536.csv")
  mouse.catSNS1536 <- ReadCatSNS1536("testdata/mouse.cat.sns.1536.csv")
  inferred.mouse.catSNS1536 <-
    TransPentanucSpectra(human.catSNS1536,
                         source.abundance = "GRCh37.genome",
                         target.abundance = "GRCm38.genome")
  expect_equal(inferred.mouse.catSNS1536, mouse.catSNS1536)
})
