context("TransformCatalogDBS144")

test_that("Transformation of a DBS 144 catalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  cat.t <- ReadCatalog("testdata/regress.cat.dbs.144.csv",
                     ref.genome = "GRCh37", region = "transcript",
                     catalog.type = "counts")

  x1 <- 
    expect_warning(TransformCatalog(cat.t, target.ref.genome = "GRCh37",
                                    target.region = "exome",
                                    target.catalog.type = "counts"))
  
  x2 <- 
    expect_warning(TransformCatalog(x1, target.ref.genome = "GRCh37",
                                    target.region = "transcript",
                                    target.catalog.type = "counts"))

  expect_equal(cat.t, x2)

  cat.density <- TransformCatalog(cat.t, target.ref.genome = "GRCh37",
                                  target.region = "transcript",
                                  target.catalog.type = "density")

  x3 <- TransformCatalog(x1, target.ref.genome = "GRCh37",
                         target.region = "exome",
                         target.catalog.type = "density")
  null.cat.density <- cat.density
  attr(null.cat.density, "region") <- NULL
  null.x3 <- x3
  attr(null.x3, "region") <- NULL
  expect_equal(null.cat.density, null.x3)
  rm(null.cat.density, null.x3)

  genome.counts.signature <-
    TransformCatalog(cat.t, target.ref.genome = "GRCh37",
                     target.region = "transcript",
                     target.catalog.type = "counts.signature")
  out <- rep(1, 4)
  expect_equal(sum(colSums(genome.counts.signature) == out), 4)

  expect_error(TransformCatalog(genome.counts.signature,
                                target.ref.genome = "GRCh37",
                                target.region = "transcript",
                                target.catalog.type = "density"))

  x4 <- TransformCatalog(genome.counts.signature,
                         target.ref.genome = "GRCh37",
                         target.region = "transcript",
                         target.catalog.type = "density.signature")

  x5 <- TransformCatalog(cat.density, 
                         target.ref.genome = "GRCh37",
                         target.region = "transcript",
                         target.catalog.type = "density.signature")
  expect_equal(x4, x5)
})
