context("TransformCatalog")

test_that("TransformCatalog function is working properly", {
cat <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")

x1 <- TransformCatalog(cat, source.abundance = "GRCh37.genome",
          target.abundance = "GRCh37.exome",
          source.type = "counts", which.n = 3)
x2 <- TransformCatalog(x1, source.abundance = "GRCh37.exome",
          target.abundance = "GRCh37.exome",
          source.type = "counts",
          target.type = "signature", which.n = 3)
x3 <- TransformCatalog(cat, source.abundance = "GRCh37.genome",
          target.abundance = "GRCh37.genome",
          source.type = "counts",
          target.type = "signature", which.n = 3)
x4 <- TransformCatalog(x3, source.abundance = "GRCh37.genome",
          target.abundance = "GRCh37.exome",
          source.type = "signature", which.n = 3)
expect_equal(x2, x4)

})
