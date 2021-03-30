context("TransformCatalogSBS192")

test_that("Transformation of a SBS 192 catalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  cat.t <- ReadCatalog("testdata/regress.cat.sbs.192.csv",
                       ref.genome = "GRCh37", 
                       region = "transcript",
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
  null.cat.density <-cat.density
  attr(null.cat.density, "region") <- NULL
  null.x3 <- x3
  attr(null.x3, "region") <- NULL
  expect_equal(null.cat.density, null.x3)
  rm(null.cat.density, null.x3)

  genome.counts.signature <-
    TransformCatalog(cat.t, target.ref.genome = "GRCh37",
                     target.region = "transcript",
                     target.catalog.type = "counts.signature")

  lapply(colSums(genome.counts.signature), FUN = function(x) {
    expect_equal(x, 1, tolerance = .Machine$double.eps^0.25)
  })

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
