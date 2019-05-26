context("TransformCatalog")

test_that("Legal transformation 1;
counts -> counts; genome counts -> exome counts,
and exome counts -> genome counts", {
  cat <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                     ref.genome = "GRCh37", region = "genome",
                     catalog.type = "counts")

  x1 <- TransformCatalog(cat, target.ref.genome = "GRCh37",
                         target.region = "exome",
                         target.catalog.type = "counts")

  x2 <- TransformCatalog(x1, target.ref.genome = "GRCh37",
                         target.region = "genome",
                         target.catalog.type = "counts")

  expect_equal(cat, x2)

})

test_that("Legal transformation 2;
counts -> density; genome counts -> density,
and genome counts -> exome count -> density", {
  cat <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                     ref.genome = "GRCh37", region = "genome",
                     catalog.type = "counts")

  x1 <- TransformCatalog(cat, target.ref.genome = "GRCh37",
                         target.region = "genome",
                         target.catalog.type = "density")

  x2 <- TransformCatalog(cat, target.ref.genome = "GRCh37",
                         target.region = "exome",
                         target.catalog.type = "counts")

  x3 <- TransformCatalog(x2, target.ref.genome = "GRCh37",
                         target.region = "exome",
                         target.catalog.type = "density")

  attr(x1, "region") <- NULL
  attr(x3, "region") <- NULL
  expect_equal(x1, x3)

})

test_that("Legal transformation 2;
counts -> density; genome GRCh37 counts -> genome GRCh37 density,
and genome GRCh37 counts -> genome GRCh38 counts -> genome GRCh38 density", {
  cat <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                     ref.genome = "GRCh37", region = "genome",
                     catalog.type = "counts")

  x1 <- TransformCatalog(cat, target.ref.genome = "GRCh37",
                         target.region = "genome",
                         target.catalog.type = "density")

  x2 <- TransformCatalog(cat, target.ref.genome = "GRCh38",
                         target.region = "genome",
                         target.catalog.type = "counts")

  x3 <- TransformCatalog(x2, target.ref.genome = "GRCh38",
                         target.region = "genome",
                         target.catalog.type = "density")

  attr(x1, "ref.genome") <- NULL
  attr(x3, "ref.genome") <- NULL
  expect_equal(x1, x3)

})

test_that("Legal transformations 1, 2, and 4;
genome GRCh37 counts -> genome GRCh38 counts,
and genome GRCh37 counts -> genome GRCh37 density -> genome GRCh38 counts", {
  cat <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                     ref.genome = "GRCh37", region = "genome",
                     catalog.type = "counts")

  x1 <- TransformCatalog(cat, target.ref.genome = "GRCh38",
                         target.region = "genome",
                         target.catalog.type = "counts")

  x2 <- TransformCatalog(cat, target.ref.genome = "GRCh37",
                         target.region = "genome",
                         target.catalog.type = "density")

  x3 <- TransformCatalog(x2, target.ref.genome = "GRCh38",
                         target.region = "genome",
                         target.catalog.type = "counts")

  expect_equal(x1, x3)

})

test_that("Legal transformation 4; density -> (genome) counts", {
  cat <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                     ref.genome = "GRCh37", region = "genome",
                     catalog.type = "counts")

  x1 <- TransformCatalog(cat, target.ref.genome = "GRCh37",
                         target.region = "genome",
                         target.catalog.type = "density")

  x2 <- TransformCatalog(x1, target.ref.genome = "GRCh37",
                         target.region = "genome",
                         target.catalog.type = "counts")

  expect_equal(cat, x2)
})

test_that("Legal transformation 5;
          density -> counts.signature, density -> density.signature", {
            cat <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                               ref.genome = "GRCh37", region = "genome",
                               catalog.type = "counts")

            x1 <- TransformCatalog(cat, target.ref.genome = "GRCh37",
                                   target.region = "genome",
                                   target.catalog.type = "density")

            x2 <- TransformCatalog(x1, target.ref.genome = "GRCh37",
                                   target.region = "genome",
                                   target.catalog.type = "counts.signature")

            x3 <- TransformCatalog(x1, target.ref.genome = "GRCh37",
                                   target.region = "genome",
                                   target.catalog.type = "density.signature")

            tmp <- apply(x1, MARGIN = 2, function(x) x/sum(x))

            expect_true(all(x3 == tmp))
            expect_true(!all(x2 == x3))
          })

test_that("Legal transformations 3 and 6;
genome counts -> exome counts.signature,
and genome counts.signature -> exome counts.signature", {
  cat <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                     ref.genome = "GRCh37", region = "genome",
                     catalog.type = "counts")

  x1 <- TransformCatalog(cat, target.ref.genome = "GRCh37",
                         target.region = "exome", target.catalog.type = "counts")

  x2 <- TransformCatalog(x1, target.ref.genome = "GRCh37",
                         target.region = "exome",
                         target.catalog.type = "counts.signature")

  x3 <- TransformCatalog(cat, target.ref.genome = "GRCh37",
                         target.region = "exome",
                         target.catalog.type = "counts.signature")
  expect_equal(x2, x3)


  x4 <- TransformCatalog(cat, target.ref.genome = "GRCh37",
                         target.region = "genome",
                         target.catalog.type = "counts.signature")

  x5 <- TransformCatalog(x4, target.ref.genome = "GRCh37",
                         target.region = "exome",
                         target.catalog.type = "counts.signature")

  expect_equal(x3, x5)

})


test_that("Legal transformation 6;
counts.signature -> denisty.signature,
and genome counts -> exome counts.signature", {
  cat <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                     ref.genome = "GRCh37", region = "genome",
                     catalog.type = "counts")

  genome.counts.signature <-
    TransformCatalog(cat, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")

  density.signature <-
    TransformCatalog(genome.counts.signature,
                     target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")

  density.signature2 <-
    TransformCatalog(cat, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")

  expect_equal(density.signature, density.signature2)

  density.cat <-
    TransformCatalog(cat, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density")

  density.signature3 <-
    TransformCatalog(density.cat, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "density.signature")

  expect_equal(density.signature, density.signature3)

  exome.counts.signature <-
    TransformCatalog(cat, target.ref.genome = "GRCh37",
                     target.region = "exome",
                     target.catalog.type = "counts.signature")

  exome.counts <-
    TransformCatalog(cat, target.ref.genome = "GRCh37",
                     target.region = "exome",
                     target.catalog.type = "counts")

  exome.counts.signature2 <-
    TransformCatalog(exome.counts, target.ref.genome = "GRCh37",
                     target.region = "exome",
                     target.catalog.type = "counts.signature")

  expect_equal(exome.counts.signature, exome.counts.signature2)

})

if (FALSE) {
  test_that("Error test: transformation of a SBS 192 catalog", {
    cat <- ReadCatalog("testdata/regress.cat.sbs.192.csv",
                       ref.genome = "GRCh37", region = "genome",
                       catalog.type = "counts")

    expect_error(TransformCatalog(cat, target.ref.genome = "GRCh37",
                                  target.region = "exome",
                                  target.catalog.type = "counts"))

    expect_error(TransformCatalog(cat, target.ref.genome = "GRCh37",
                                  target.region = "genome",
                                  target.catalog.type = "density"))

    genome.counts.signature <-
      TransformCatalog(cat, target.ref.genome = "GRCh37",
                       target.region = "genome",
                       target.catalog.type = "counts.signature")
    out <- rep(1, 4)
    expect_equal(sum(colSums(genome.counts.signature) == out), 4)

    expect_error(TransformCatalog(genome.counts.signature,
                                  target.ref.genome = "GRCh37",
                                  target.region = "genome",
                                  target.catalog.type = "density"))

    expect_error(TransformCatalog(genome.counts.signature,
                                  target.ref.genome = "GRCh37",
                                  target.region = "genome",
                                  target.catalog.type = "density.signature"))
  })
}


if (FALSE) {
  test_that("Transformation of a DBS 144 catalog", {
    cat <- ReadCatalog("testdata/regress.cat.dbs.144.csv",
                       ref.genome = "GRCh37", region = "transcript",
                       catalog.type = "counts")

    cat.exome <- TransformCatalog(cat, target.ref.genome = "GRCh37",
                                  target.region = "exome",
                                  target.catalog.type = "counts")

    expect_error(TransformCatalog(cat, target.ref.genome = "GRCh37",
                                  target.region = "genome",
                                  target.catalog.type = "density"))

    genome.counts.signature <-
      TransformCatalog(cat, target.ref.genome = "GRCh37",
                       target.region = "genome",
                       target.catalog.type = "counts.signature")
    out <- rep(1, 4)
    expect_equal(sum(colSums(genome.counts.signature) == out), 4)

    expect_error(TransformCatalog(genome.counts.signature,
                                  target.ref.genome = "GRCh37",
                                  target.region = "genome",
                                  target.catalog.type = "density"))

    expect_error(TransformCatalog(genome.counts.signature,
                                  target.ref.genome = "GRCh37",
                                  target.region = "genome",
                                  target.catalog.type = "density.signature"))
  })
}


test_that("going from density to density,
          density.signature to density.signature", {
            cat <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                               ref.genome = "GRCh37", region = "genome",
                               catalog.type = "counts")

            x1 <- TransformCatalog(cat, target.ref.genome = "GRCh37",
                                   target.region = "genome",
                                   target.catalog.type = "density")


            x2 <- TransformCatalog(x1, target.ref.genome = "GRCh37",
                                   target.region = "genome",
                                   target.catalog.type = "density")

            x3 <- TransformCatalog(x1, target.ref.genome = "GRCh37",
                                   target.region = "genome",
                                   target.catalog.type = "density.signature")

            x4 <- TransformCatalog(x3, target.ref.genome = "GRCh37",
                                   target.region = "genome",
                                   target.catalog.type = "density.signature")

            expect_equal(x1, x2)
            expect_equal(x3, x4)
          })

test_that("Error test: counts.singature -> counts or density,
            error message expected", {
              cat <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                                 ref.genome = "GRCh37", region = "genome",
                                 catalog.type = "counts")

              x1 <- TransformCatalog(cat, target.ref.genome = "GRCh37",
                                     target.region = "exome",
                                     target.catalog.type = "counts.signature")

              expect_error(TransformCatalog(x1, target.ref.genome = "GRCh37",
                                            target.region = "exome",
                                            target.catalog.type = "counts"))

              expect_error(TransformCatalog(x1, target.ref.genome = "GRCh37",
                                            target.region = "exome",
                                            target.catalog.type = "density"))
            })

