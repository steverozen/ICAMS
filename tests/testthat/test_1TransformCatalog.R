context("TransformCatalog")

test_that("Legal transformation 1;
counts -> counts; genome counts -> exome counts,
and exome counts -> genome counts", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  cat <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                     ref.genome = "GRCh37", region = "genome",
                     catalog.type = "counts")

  x1 <- 
    expect_warning(TransformCatalog(cat, target.ref.genome = "GRCh37",
                         target.region = "exome",
                         target.catalog.type = "counts"))
  

  x2 <- 
    expect_warning(TransformCatalog(x1, target.ref.genome = "GRCh37",
                         target.region = "genome",
                         target.catalog.type = "counts"))

  expect_equal(cat, x2)

})

test_that("Legal transformation 2;
counts -> density; genome counts -> density,
and genome counts -> exome count -> density", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  cat <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                     ref.genome = "GRCh37", region = "genome",
                     catalog.type = "counts")

  x1 <- TransformCatalog(cat, target.ref.genome = "GRCh37",
                         target.region = "genome",
                         target.catalog.type = "density")

  x2 <- 
    expect_warning(TransformCatalog(cat, target.ref.genome = "GRCh37",
                                    target.region = "exome",
                                    target.catalog.type = "counts"))

  x3 <- TransformCatalog(x2, target.ref.genome = "GRCh37",
                         target.region = "exome",
                         target.catalog.type = "density")

  attr(x1, "region") <- NULL
  attr(x3, "region") <- NULL
  expect_equal(x1, x3)

})

test_that("Legal transformation 3;
counts -> density; genome GRCh37 counts -> genome GRCh37 density,
and genome GRCh37 counts -> genome GRCh38 counts -> genome GRCh38 density", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.UCSC.hg38"))
  cat <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                     ref.genome = "GRCh37", region = "genome",
                     catalog.type = "counts")

  x1 <- TransformCatalog(cat, target.ref.genome = "GRCh37",
                         target.region = "genome",
                         target.catalog.type = "density")

  x2 <- 
    expect_warning(TransformCatalog(cat, target.ref.genome = "GRCh38",
                                    target.region = "genome",
                                    target.catalog.type = "counts"))

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
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.UCSC.hg38"))
  cat <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                     ref.genome = "GRCh37", region = "genome",
                     catalog.type = "counts")

  x1 <- 
    expect_warning(TransformCatalog(cat, target.ref.genome = "GRCh38",
                                    target.region = "genome",
                                    target.catalog.type = "counts"))
  
  x2 <- TransformCatalog(cat, target.ref.genome = "GRCh37",
                         target.region = "genome",
                         target.catalog.type = "density")
  
  x3 <- TransformCatalog(x2, target.ref.genome = "GRCh38",
                         target.region = "genome",
                         target.catalog.type = "counts")

  expect_equal(x1, x3)

})

test_that("Legal transformation 4; density -> (genome) counts", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
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
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
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
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  cat <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                     ref.genome = "GRCh37", region = "genome",
                     catalog.type = "counts")

  x1 <- 
    expect_warning(TransformCatalog(cat, target.ref.genome = "GRCh37",
                                    target.region = "exome", 
                                    target.catalog.type = "counts"))

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
counts.signature -> density.signature,
and genome counts -> exome counts.signature", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
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
    expect_warning(TransformCatalog(cat, target.ref.genome = "GRCh37",
                                    target.region = "exome",
                                    target.catalog.type = "counts"))

  exome.counts.signature2 <-
    TransformCatalog(exome.counts, target.ref.genome = "GRCh37",
                     target.region = "exome",
                     target.catalog.type = "counts.signature")

  expect_equal(exome.counts.signature, exome.counts.signature2)

})

test_that("going from density to density,
          density.signature to density.signature", {
            skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
            stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
            cat <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                               ref.genome = "GRCh37", region = "genome",
                               catalog.type = "counts")

            x1 <- TransformCatalog(cat, target.ref.genome = "GRCh37",
                                   target.region = "genome",
                                   target.catalog.type = "density")


            expect_warning(
              x2 <- TransformCatalog(x1, target.ref.genome = "GRCh37",
                                   target.region = "genome",
                                   target.catalog.type = "density"))

            x3 <- TransformCatalog(x1, target.ref.genome = "GRCh37",
                                   target.region = "genome",
                                   target.catalog.type = "density.signature")

            expect_warning(x4 <- TransformCatalog(x3, target.ref.genome = "GRCh37",
                                   target.region = "genome",
                                   target.catalog.type = "density.signature"))

            expect_equal(x1, x2)
            expect_equal(x3, x4)
          })

test_that("Error test: counts.singature -> counts or density,
            error message expected", {
              skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
              stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
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

test_that("Transform a counts catalog with NULL abundance,
          to counts.signature catalog", {
            skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
            skip_if("" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38"))
            stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
            stopifnot(requireNamespace("BSgenome.Hsapiens.UCSC.hg38"))
            catSBS96.counts <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                                           ref.genome = "hg19",
                                           catalog.type = "counts")
            catSBS96.counts.sig <- 
              TransformCatalog(catSBS96.counts,
                               target.catalog.type = "counts.signature")
            expect_equal(sum(colSums(catSBS96.counts.sig) == rep(1, 4)), 4)
            
            expect_error(
              TransformCatalog(catSBS96.counts, target.ref.genome = "hg38",
                               target.catalog.type = "counts.signature")
            )
            
            
            catSBS192.counts <- ReadCatalog("testdata/regress.cat.sbs.192.csv",
                                            catalog.type = "counts")
            catSBS192.counts.sig <- 
              TransformCatalog(catSBS192.counts,
                               target.catalog.type = "counts.signature")
            expect_equal(sum(colSums(catSBS192.counts.sig) == rep(1, 4)), 4)
            expect_error(
              TransformCatalog(catSBS192.counts, target.ref.genome = "hg38",
                               target.catalog.type = "counts.signature")
            )
            
            
            catSBS1536.counts <- ReadCatalog("testdata/regress.cat.sbs.1536.csv",
                                             catalog.type = "counts")
            catSBS1536.counts.sig <- 
              TransformCatalog(catSBS1536.counts,
                               target.catalog.type = "counts.signature")
            expect_equal(sum(colSums(catSBS1536.counts.sig) == rep(1, 4)), 4)
            expect_error(
              TransformCatalog(catSBS1536.counts, target.ref.genome = "hg38",
                               target.catalog.type = "counts.signature")
            )
            
            
            catDBS78.counts <- ReadCatalog("testdata/regress.cat.dbs.78.csv",
                                           catalog.type = "counts")
            catDBS78.counts.sig <- 
              TransformCatalog(catDBS78.counts,
                               target.catalog.type = "counts.signature")
            expect_equal(sum(colSums(catDBS78.counts.sig) == rep(1, 4)), 4)
            expect_error(
              TransformCatalog(catDBS78.counts, target.ref.genome = "hg38",
                               target.catalog.type = "counts.signature")
            )
            
            catDBS144.counts <- ReadCatalog("testdata/regress.cat.dbs.144.csv",
                                            catalog.type = "counts")
            catDBS144.counts.sig <- 
              TransformCatalog(catDBS144.counts,
                               target.catalog.type = "counts.signature")
            expect_equal(sum(colSums(catDBS144.counts.sig) == rep(1, 4)), 4)
            expect_error(
              TransformCatalog(catDBS144.counts, target.ref.genome = "hg38",
                               target.catalog.type = "counts.signature")
            )
            
            
            catDBS136.counts <- ReadCatalog("testdata/regress.cat.dbs.136.csv",
                                            catalog.type = "counts")
            catDBS136.counts.sig <- 
              TransformCatalog(catDBS136.counts,
                               target.catalog.type = "counts.signature")
            expect_equal(sum(colSums(catDBS136.counts.sig) == rep(1, 4)), 4)
            expect_error(
              TransformCatalog(catDBS136.counts, target.ref.genome = "hg38",
                               target.catalog.type = "counts.signature")
            )
            
            
            catID.counts <- ReadCatalog("testdata/BTSG_WGS_PCAWG.indels.csv",
                                        catalog.type = "counts")
            catID.counts.sig <- 
              TransformCatalog(catID.counts,
                               target.catalog.type = "counts.signature")
            expect_equal(sum(colSums(catID.counts.sig) == rep(1, 35)), 35)
            expect_error(
              TransformCatalog(catID.counts, target.ref.genome = "hg38",
                               target.catalog.type = "counts.signature")
            )
})

