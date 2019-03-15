context("TransformCatalog")

test_that("TransformCatalog genome counts -> exome signature,
          and genome signature -> exome signature", {
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


test_that("TransformCatalog genome-count signature -> denisty signature,
          and genome counts -> exome-count signature", {
            cat <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")

            genome.count.signature <-
              TransformCatalog(cat, source.abundance = "GRCh37.genome",
                               target.abundance = "GRCh37.genome",
                               source.type = "counts",
                               target.type = "signature", which.n = 3)
            density.signature <-
              TransformCatalog(genome.count.signature,
                               source.abundance = "GRCh37.genome",
                               target.abundance = NULL,
                               source.type = "signature",
                               target.type = "signature",
                               which.n = 3)
            density.cat <-
              TransformCatalog(cat, source.abundance = "GRCh37.genome",
                               target.abundance = NULL,
                               source.type = "counts",
                               target.type = "density",
                               which.n = 3)

            density.signature2 <-
              TransformCatalog(density.cat, source.abundance = NULL,
                               source.type = "density",
                               target.type = "signature",
                               which.n = 3)

            exome.count.signature <-
              TransformCatalog(cat, source.abundance = "GRCh37.genome",
                               target.abundance = "GRCh37.exome",
                               source.type = "counts",
                               target.type = "signature", which.n = 3)

            exome.counts <-
              TransformCatalog(cat, source.abundance = "GRCh37.genome",
                               target.abundance = "GRCh37.exome",
                               source.type = "counts",
                               target.type = "counts", which.n = 3)

            exome.count.signature2 <-
              TransformCatalog(exome.counts, source.abundance = "GRCh37.exome",
                               target.abundance = "GRCh37.exome",
                               source.type = "counts",
                               target.type = "signature", which.n = 3)

            expect_equal(density.signature2, density.signature)

            })

test_that("TransformCatalog genome counts -> exome counts,
          and exome counts -> genome counts", {
            cat <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")

            x1 <- TransformCatalog(cat, source.abundance = "GRCh37.genome",
                                   target.abundance = "GRCh37.exome",
                                   source.type = "counts", which.n = 3)
            x2 <- TransformCatalog(x1, source.abundance = "GRCh37.exome",
                                   target.abundance = "GRCh37.genome",
                                   source.type = "counts", which.n = 3)
            expect_equal(cat, x2)

          })

test_that("TransformCatalog genome counts -> density,
           and genome counts -> exome count ->density", {
            cat <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")

            expect_error(TransformCatalog(cat, source.abundance = "GRCh37.genome",
                                          target.abundance = "GRCh37.genome",
                                          source.type = "counts",
                                          target.type = "density", which.n = 3))

            x1 <- TransformCatalog(cat, source.abundance = "GRCh37.genome",
                                   target.abundance = NULL,
                                   source.type = "counts",
                                   target.type = "density", which.n = 3)

            x2 <- TransformCatalog(cat, source.abundance = "GRCh37.genome",
                                   target.abundance = "GRCh37.exome",
                                   source.type = "counts", which.n = 3)

            x3 <- TransformCatalog(x2, source.abundance = "GRCh37.exome",
                                   target.abundance = NULL,
                                   source.type = "counts",
                                   target.type = "density", which.n = 3)
            expect_equal(x1, x3)

          })

test_that("TransformCatalog genome GRCh37 counts -> genome GRCh37 density,
           and genome GRCh37 counts -> genome GRCh38 count -> genome GRCh38 density", {
             cat <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")

             x1 <- TransformCatalog(cat, source.abundance = "GRCh37.genome",
                                    target.abundance = NULL,
                                    source.type = "counts",
                                    target.type = "density", which.n = 3)

             x2 <- TransformCatalog(cat, source.abundance = "GRCh37.genome",
                                    target.abundance = "GRCh38.genome",
                                    source.type = "counts",
                                    target.type = "counts", which.n = 3)

             x3 <- TransformCatalog(x2, source.abundance = "GRCh38.genome",
                                    target.abundance = NULL,
                                    source.type = "counts",
                                    target.type = "density", which.n = 3)

             expect_equal(x1, x3)

           })

test_that("TransformCatalog genome GRCh37 counts -> genome GRCh38 counts,
           and genome GRCh37 counts -> genome GRCh37 density -> genome GRCh38 counts", {
             cat <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")

             x1 <- TransformCatalog(cat, source.abundance = "GRCh37.genome",
                                    target.abundance = "GRCh38.genome",
                                    source.type = "counts",
                                    target.type = "counts", which.n = 3)

             x2 <- TransformCatalog(cat, source.abundance = "GRCh37.genome",
                                    target.abundance = NULL,
                                    source.type = "counts",
                                    target.type = "density", which.n = 3)

             x3 <- TransformCatalog(x2, source.abundance = NULL,
                                    target.abundance = "GRCh38.genome",
                                    source.type = "density",
                                    target.type = "counts", which.n = 3)

             expect_equal(x1, x3)

           })

test_that("TransformCatalog function: going from density to genome counts", {
  cat <- ReadCatSNS96("testdata/regress.cat.sns.96.csv")

  x1 <- TransformCatalog(cat, source.abundance = "GRCh37.genome",
                         source.type = "counts",
                         target.type = "density", which.n = 3)

  x2 <- TransformCatalog(x1, source.abundance = NULL,
                         source.type = "density",
                         target.abundance = "GRCh37.genome",
                         target.type = "counts", which.n = 3)

  expect_equal(cat, x2)
})

test_that("TransformCatalog function: transformation of a SNS 192 catalog", {
            cat <- ReadCatSNS192("testdata/regress.cat.sns.192.csv")

            expect_error(TransformCatalog(cat, source.abundance = "GRCh37.genome",
                                          target.abundance = "GRCh37.exome",
                                          source.type = "counts",
                                          target.type = "counts", which.n = 3))

            expect_error(TransformCatalog(cat, source.abundance = "GRCh37.genome",
                                          source.type = "counts",
                                          target.type = "density", which.n = 3))

            genome.count.signature <-
              TransformCatalog(cat, source.abundance = "GRCh37.genome",
                               target.abundance = "GRCh37.genome",
                               source.type = "counts",
                               target.type = "signature", which.n = 3)
            out <- rep(1, 4)
            expect_equal(sum(colSums(genome.count.signature) == out), 4)

            expect_error(TransformCatalog(genome.count.signature,
                                          source.abundance = "GRCh37.genome",
                                          target.abundance = NULL,
                                          source.type = "signature",
                                          target.type = "density",
                                          which.n = 3))

            expect_error(TransformCatalog(genome.count.signature,
                                          source.abundance = "GRCh37.genome",
                                          target.abundance = NULL,
                                          source.type = "signature",
                                          target.type = "signature",
                                          which.n = 3))
          })

test_that("TransformCatalog function: transformation of a DNS 144 catalog", {
  cat <- ReadCatDNS144("testdata/regress.cat.dns.144.csv")

  expect_error(TransformCatalog(cat, source.abundance = "GRCh37.genome",
                                target.abundance = "GRCh37.exome",
                                source.type = "counts",
                                target.type = "counts", which.n = 2))

  expect_error(TransformCatalog(cat, source.abundance = "GRCh37.genome",
                                source.type = "counts",
                                target.type = "density", which.n = 2))

  genome.count.signature <-
    TransformCatalog(cat, source.abundance = "GRCh37.genome",
                     target.abundance = "GRCh37.genome",
                     source.type = "counts",
                     target.type = "signature", which.n = 2)
  out <- rep(1, 4)
  expect_equal(sum(colSums(genome.count.signature) == out), 4)

  expect_error(TransformCatalog(genome.count.signature,
                                source.abundance = "GRCh37.genome",
                                target.abundance = NULL,
                                source.type = "signature",
                                target.type = "density",
                                which.n = 2))

  expect_error(TransformCatalog(genome.count.signature,
                                source.abundance = "GRCh37.genome",
                                target.abundance = NULL,
                                source.type = "signature",
                                target.type = "signature",
                                which.n = 2))
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
