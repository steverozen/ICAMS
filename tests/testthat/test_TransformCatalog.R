context("TransformCatalog")

test_that("TransformCatalog function: going from genome counts to exome signature,
          from genome signature to exome signature", {
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


test_that("TransformCatalog function: going from genome counts to exome counts,
          from exome counts back to genome counts", {
            cat <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")

            x1 <- TransformCatalog(cat, source.abundance = "GRCh37.genome",
                                   target.abundance = "GRCh37.exome",
                                   source.type = "counts", which.n = 3)
            x2 <- TransformCatalog(x1, source.abundance = "GRCh37.exome",
                                   target.abundance = "GRCh37.genome",
                                   source.type = "counts", which.n = 3)
            expect_equal(cat, x2)

          })

test_that("TransformCatalog function: going from genome counts to genome density,
          from genome counts to exome density", {
            cat <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")

            x1 <- TransformCatalog(cat, source.abundance = "GRCh37.genome",
                                   target.abundance = "GRCh37.genome",
                                   source.type = "counts",
                                   target.type = "density", which.n = 3)
            x2 <- TransformCatalog(cat, source.abundance = "GRCh37.genome",
                                   target.abundance = "GRCh37.exome",
                                   source.type = "counts", which.n = 3)
            x3 <- TransformCatalog(x2, source.abundance = "GRCh37.exome",
                                   target.abundance = "GRCh37.exome",
                                   source.type = "counts",
                                   target.type = "density", which.n = 3)
            expect_equal(x1, x3)

          })

test_that("TransformCatalog function: going from density to density,
          error message expected", {
            cat <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")

            x1 <- TransformCatalog(cat, source.abundance = "GRCh37.genome",
                                   source.type = "counts",
                                   target.type = "density", which.n = 3)

            expect_error(TransformCatalog(x1, source.abundance = "GRCh37.genome",
                                          source.type = "density",
                                          target.type = "density", which.n = 3))
          })

test_that("TransformCatalog function: going from singature to counts or density,
          error message expected", {
            cat <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")

            x1 <- TransformCatalog(cat, source.abundance = "GRCh37.genome",
                                   target.abundance = "GRCh37.exome",
                                   source.type = "counts",
                                   target.type = "signature", which.n = 3)

            expect_error(TransformCatalog(x1, source.abundance = "GRCh37.exome",
                                          target.abundance = "GRCh37.genome",
                                          source.type = "signature",
                                          target.type = "counts", which.n = 3))

            expect_error(TransformCatalog(x1, source.abundance = "GRCh37.exome",
                                          target.abundance = "GRCh37.genome",
                                          source.type = "signature",
                                          target.type = "density", which.n = 3))
          })

test_that("TransformCatalog function: specifying the wrong which.n argument,
          error message expected", {
            cat <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")

            expect_error(TransformCatalog(cat, source.abundance = "GRCh37.genome",
                                          target.abundance = "GRCh37.exome",
                                          source.type = "counts",
                                          target.type = "signature", which.n = 4))
          })

test_that("TransformCatalog function: specifying the wrong abundance,
          error message expected", {
            cat <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")

            expect_error(TransformCatalog(cat, source.abundance = "GRC37.genome",
                                          target.abundance = "GRCh37.exome",
                                          source.type = "counts",
                                          target.type = "signature", which.n = 3))
          })

test_that("TransformCatalog function: abundance is not specified,
          error message expected", {
            cat <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")

            expect_error(TransformCatalog(cat, source.abundance = "GRC37.genome",
                                          source.type = "counts",
                                          target.type = "signature", which.n = 3),
                         "Please specify the target.abundance", fixed = TRUE)
          })
