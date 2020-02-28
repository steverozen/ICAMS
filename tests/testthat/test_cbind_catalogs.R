context("cbind methods for catalogs")

test_that("cbind method for SBS96Catalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  catalog.counts <- ReadCatalog("testdata/regress.cat.sbs.96.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)

  cat1 <- catalog.counts[, 1, drop = FALSE]
  cat2 <- catalog.counts[, 2, drop = FALSE]
  cat3 <- cbind(cat1, cat2)
  out <- PlotCatalogToPdf(cat3, file = file.path(tempdir(), "test0.pdf"))
  expect_equal(out, TRUE)
  
  # Test when performing cbind operation to catalogs with different "ref.genome"
  # attribute
  cat4 <- catalog.counts[, 3, drop = FALSE]
  
  # When some catalog has NULL "ref.genome" attribute
  attr(cat4, "ref.genome") <- NULL
  combined.cat <- cbind(cat1, cat4)
  expect_equal(attr(combined.cat, "ref.genome"), "mixed")
  out1 <- 
    PlotCatalogToPdf(combined.cat, file = file.path(tempdir(), "test1.pdf"))
  expect_equal(out1, TRUE)
  expect_error(TransformCatalog(combined.cat, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "ref.genome" attribute
  cat5 <- cat4
  combined.cat1 <- cbind(cat4, cat5)
  expect_equal(attr(combined.cat1, "ref.genome"), NULL)
  out2 <- 
    PlotCatalogToPdf(combined.cat1, file = file.path(tempdir(), "test2.pdf"))
  expect_equal(out2, TRUE)
  t.cat1 <- TransformCatalog(combined.cat1, target.ref.genome = "hg19",
                            target.catalog.type = "counts.signature")
  
  # When catalogs have non-NULL different "ref.genome" attribute
  cat6 <- cat1
  if (!"" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38")) {
    attr(cat6, "ref.genome") <-
      BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    combined.cat2 <- cbind(cat2, cat6)
    expect_equal(attr(combined.cat2, "ref.genome"), "mixed")
    out3 <- 
      PlotCatalogToPdf(combined.cat2, file = file.path(tempdir(), "test3.pdf"))
    expect_equal(out3, TRUE)
    unlink(file.path(tempdir(), "test3.pdf"))
    expect_error(TransformCatalog(combined.cat2, target.ref.genome = "hg19",
                              target.catalog.type = "counts.signature"))
  }
  
  # Test when performing cbind operation to catalogs with different
  # "catalog.type" attribute
  cat7 <- cat8 <- TransformCatalog(cat1, target.catalog.type = "density")
  
  # When catalogs have non-NULL different "catalog.type" attribute
  expect_error(cbind(cat2, cat7))
  
  # When some catalog has NULL "catalog.type" attribute
  attr(cat7, "catalog.type") <- NULL
  expect_error(cbind(cat2, cat7))
  
  # When all catalogs have NULL "catalog.type" attribute
  attr(cat8, "catalog.type") <- NULL
  expect_error(cbind(cat7, cat8))
  expect_error(TransformCatalog(cat8, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # Test when performing cbind operation to catalogs with different "abundance"
  # attribute
  cat9 <- cat10 <- cat1
  
  # When catalogs have non-NULL different "abundance" attribute
  attr(cat9, "abundance") <- abundance.3bp.genome.unstranded.GRCh38
  combined.cat4 <- cbind(cat1, cat9)
  expect_equal(attr(combined.cat4, "abundance"), "mixed")
  out4 <- 
    PlotCatalogToPdf(combined.cat4, file = file.path(tempdir(), "test4.pdf"))
  expect_equal(out4, TRUE)
  expect_error(TransformCatalog(combined.cat4, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When some catalog have NULL "abundance" attribute
  attr(cat9, "abundance") <- NULL
  combined.cat5 <- cbind(cat1, cat9)
  expect_equal(attr(combined.cat5, "abundance"), "mixed")
  out5 <- 
    PlotCatalogToPdf(combined.cat5, file = file.path(tempdir(), "test5.pdf"))
  expect_equal(out5, TRUE)
  expect_error(TransformCatalog(combined.cat5, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "abundance" attribute
  attr(cat10, "abundance") <- NULL
  combined.cat6 <- cbind(cat9, cat10)
  expect_equal(attr(combined.cat6, "abundance"), NULL)
  out6 <- 
    PlotCatalogToPdf(combined.cat6, file = file.path(tempdir(), "test6.pdf"))
  expect_equal(out6, TRUE)
  expect_error(TransformCatalog(combined.cat6, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # Test when performing cbind operation to catalogs with different "region"
  # attribute
  cat11 <- cat12 <- cat1
  
  # When catalogs have non-NULL different "region" attribute
  attr(cat11, "region") <- "exome"
  combined.cat7 <- cbind(cat1, cat11)
  expect_equal(attr(combined.cat7, "region"), "mixed")
  out7 <- 
    PlotCatalogToPdf(combined.cat7, file = file.path(tempdir(), "test7.pdf"))
  expect_equal(out7, TRUE)
  expect_error(TransformCatalog(combined.cat7, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When some catalog have NULL "region" attribute
  attr(cat11, "region") <- NULL
  combined.cat8 <- cbind(cat1, cat11)
  expect_equal(attr(combined.cat8, "region"), "mixed")
  out8 <- 
    PlotCatalogToPdf(combined.cat8, file = file.path(tempdir(), "test8.pdf"))
  expect_equal(out8, TRUE)
  expect_error(TransformCatalog(combined.cat8, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "region" attribute
  attr(cat12, "region") <- NULL
  combined.cat9 <- cbind(cat11, cat12)
  expect_equal(attr(combined.cat9, "region"), NULL)
  out9 <- 
    PlotCatalogToPdf(combined.cat9, file = file.path(tempdir(), "test9.pdf"))
  expect_equal(out9, TRUE)
  expect_error(TransformCatalog(combined.cat9, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
    
  graphics.off()
  for (i in 0:9) {
    path <- file.path(tempdir(), paste0("test", i, ".pdf"))
    unlink(path)
  }
})

test_that("cbind method for SBS192Catalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  catalog.counts <- ReadCatalog("testdata/regress.cat.sbs.192.csv",
                                ref.genome = "GRCh37",
                                region = "transcript",
                                catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)

  cat1 <- catalog.counts[, 1, drop = FALSE]
  cat2 <- catalog.counts[, 2, drop = FALSE]
  cat3 <- cbind(cat1, cat2)
  out <- PlotCatalogToPdf(cat3, file = file.path(tempdir(), "test0.pdf"))
  expect_equal(out, TRUE)
  
  # Test when performing cbind operation to catalogs with different "ref.genome"
  # attribute
  cat4 <- catalog.counts[, 3, drop = FALSE]
  
  # When some catalog has NULL "ref.genome" attribute
  attr(cat4, "ref.genome") <- NULL
  combined.cat <- cbind(cat1, cat4)
  expect_equal(attr(combined.cat, "ref.genome"), "mixed")
  out1 <- 
    PlotCatalogToPdf(combined.cat, file = file.path(tempdir(), "test1.pdf"))
  expect_equal(out1, TRUE)
  expect_error(TransformCatalog(combined.cat, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "ref.genome" attribute
  cat5 <- cat4
  combined.cat1 <- cbind(cat4, cat5)
  expect_equal(attr(combined.cat1, "ref.genome"), NULL)
  out2 <- 
    PlotCatalogToPdf(combined.cat1, file = file.path(tempdir(), "test2.pdf"))
  expect_equal(out2, TRUE)
  t.cat1 <- TransformCatalog(combined.cat1, target.ref.genome = "hg19",
                             target.catalog.type = "counts.signature")
  
  # When catalogs have non-NULL different "ref.genome" attribute
  cat6 <- cat1
  if (!"" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38")) {
    attr(cat6, "ref.genome") <-
      BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    combined.cat2 <- cbind(cat2, cat6)
    expect_equal(attr(combined.cat2, "ref.genome"), "mixed")
    out3 <- 
      PlotCatalogToPdf(combined.cat2, file = file.path(tempdir(), "test3.pdf"))
    expect_equal(out3, TRUE)
    unlink(file.path(tempdir(), "test3.pdf"))
    expect_error(TransformCatalog(combined.cat2, target.ref.genome = "hg19",
                                  target.catalog.type = "counts.signature"))
  }
  
  # Test when performing cbind operation to catalogs with different
  # "catalog.type" attribute
  cat7 <- cat8 <- TransformCatalog(cat1, target.catalog.type = "density")
  
  # When catalogs have non-NULL different "catalog.type" attribute
  expect_error(cbind(cat2, cat7))
  
  # When some catalog has NULL "catalog.type" attribute
  attr(cat7, "catalog.type") <- NULL
  expect_error(cbind(cat2, cat7))
  
  # When all catalogs have NULL "catalog.type" attribute
  attr(cat8, "catalog.type") <- NULL
  expect_error(cbind(cat7, cat8))
  expect_error(TransformCatalog(cat8, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # Test when performing cbind operation to catalogs with different "abundance"
  # attribute
  cat9 <- cat10 <- cat1
  
  # When catalogs have non-NULL different "abundance" attribute
  attr(cat9, "abundance") <- abundance.3bp.genome.unstranded.GRCh38
  combined.cat4 <- cbind(cat1, cat9)
  expect_equal(attr(combined.cat4, "abundance"), "mixed")
  out4 <- 
    PlotCatalogToPdf(combined.cat4, file = file.path(tempdir(), "test4.pdf"))
  expect_equal(out4, TRUE)
  expect_error(TransformCatalog(combined.cat4, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When some catalog have NULL "abundance" attribute
  attr(cat9, "abundance") <- NULL
  combined.cat5 <- cbind(cat1, cat9)
  expect_equal(attr(combined.cat5, "abundance"), "mixed")
  out5 <- 
    PlotCatalogToPdf(combined.cat5, file = file.path(tempdir(), "test5.pdf"))
  expect_equal(out5, TRUE)
  expect_error(TransformCatalog(combined.cat5, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "abundance" attribute
  attr(cat10, "abundance") <- NULL
  combined.cat6 <- cbind(cat9, cat10)
  expect_equal(attr(combined.cat6, "abundance"), NULL)
  out6 <- 
    PlotCatalogToPdf(combined.cat6, file = file.path(tempdir(), "test6.pdf"))
  expect_equal(out6, TRUE)
  expect_error(TransformCatalog(combined.cat6, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # Test when performing cbind operation to catalogs with different "region"
  # attribute
  cat11 <- cat12 <- cat1
  
  # When catalogs have non-NULL different "region" attribute
  attr(cat11, "region") <- "exome"
  combined.cat7 <- cbind(cat1, cat11)
  expect_equal(attr(combined.cat7, "region"), "mixed")
  out7 <- 
    PlotCatalogToPdf(combined.cat7, file = file.path(tempdir(), "test7.pdf"))
  expect_equal(out7, TRUE)
  expect_error(TransformCatalog(combined.cat7, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When some catalog have NULL "region" attribute
  attr(cat11, "region") <- NULL
  combined.cat8 <- cbind(cat1, cat11)
  expect_equal(attr(combined.cat8, "region"), "mixed")
  out8 <- 
    PlotCatalogToPdf(combined.cat8, file = file.path(tempdir(), "test8.pdf"))
  expect_equal(out8, TRUE)
  expect_error(TransformCatalog(combined.cat8, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "region" attribute
  attr(cat12, "region") <- NULL
  combined.cat9 <- cbind(cat11, cat12)
  expect_equal(attr(combined.cat9, "region"), NULL)
  out9 <- 
    PlotCatalogToPdf(combined.cat9, file = file.path(tempdir(), "test9.pdf"))
  expect_equal(out9, TRUE)
  expect_error(TransformCatalog(combined.cat9, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  graphics.off()
  for (i in 0:9) {
    path <- file.path(tempdir(), paste0("test", i, ".pdf"))
    unlink(path)
  }
})

test_that("cbind method for SBS1536Catalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  catalog.counts <- ReadCatalog("testdata/regress.cat.sbs.1536.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)

  cat1 <- catalog.counts[, 1, drop = FALSE]
  cat2 <- catalog.counts[, 2, drop = FALSE]
  cat3 <- cbind(cat1, cat2)
  out <- PlotCatalogToPdf(cat3, file = file.path(tempdir(), "test0.pdf"))
  expect_equal(out, TRUE)
  
  # Test when performing cbind operation to catalogs with different "ref.genome"
  # attribute
  cat4 <- catalog.counts[, 3, drop = FALSE]
  
  # When some catalog has NULL "ref.genome" attribute
  attr(cat4, "ref.genome") <- NULL
  combined.cat <- cbind(cat1, cat4)
  expect_equal(attr(combined.cat, "ref.genome"), "mixed")
  out1 <- 
    PlotCatalogToPdf(combined.cat, file = file.path(tempdir(), "test1.pdf"))
  expect_equal(out1, TRUE)
  expect_error(TransformCatalog(combined.cat, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "ref.genome" attribute
  cat5 <- cat4
  combined.cat1 <- cbind(cat4, cat5)
  expect_equal(attr(combined.cat1, "ref.genome"), NULL)
  out2 <- 
    PlotCatalogToPdf(combined.cat1, file = file.path(tempdir(), "test2.pdf"))
  expect_equal(out2, TRUE)
  t.cat1 <- TransformCatalog(combined.cat1, target.ref.genome = "hg19",
                             target.catalog.type = "counts.signature")
  
  # When catalogs have non-NULL different "ref.genome" attribute
  cat6 <- cat1
  if (!"" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38")) {
    attr(cat6, "ref.genome") <-
      BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    combined.cat2 <- cbind(cat2, cat6)
    expect_equal(attr(combined.cat2, "ref.genome"), "mixed")
    out3 <- 
      PlotCatalogToPdf(combined.cat2, file = file.path(tempdir(), "test3.pdf"))
    expect_equal(out3, TRUE)
    unlink(file.path(tempdir(), "test3.pdf"))
    expect_error(TransformCatalog(combined.cat2, target.ref.genome = "hg19",
                                  target.catalog.type = "counts.signature"))
  }
  
  # Test when performing cbind operation to catalogs with different
  # "catalog.type" attribute
  cat7 <- cat8 <- TransformCatalog(cat1, target.catalog.type = "density")
  
  # When catalogs have non-NULL different "catalog.type" attribute
  expect_error(cbind(cat2, cat7))
  
  # When some catalog has NULL "catalog.type" attribute
  attr(cat7, "catalog.type") <- NULL
  expect_error(cbind(cat2, cat7))
  
  # When all catalogs have NULL "catalog.type" attribute
  attr(cat8, "catalog.type") <- NULL
  expect_error(cbind(cat7, cat8))
  expect_error(TransformCatalog(cat8, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # Test when performing cbind operation to catalogs with different "abundance"
  # attribute
  cat9 <- cat10 <- cat1
  
  # When catalogs have non-NULL different "abundance" attribute
  attr(cat9, "abundance") <- abundance.3bp.genome.unstranded.GRCh38
  combined.cat4 <- cbind(cat1, cat9)
  expect_equal(attr(combined.cat4, "abundance"), "mixed")
  out4 <- 
    PlotCatalogToPdf(combined.cat4, file = file.path(tempdir(), "test4.pdf"))
  expect_equal(out4, TRUE)
  expect_error(TransformCatalog(combined.cat4, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When some catalog have NULL "abundance" attribute
  attr(cat9, "abundance") <- NULL
  combined.cat5 <- cbind(cat1, cat9)
  expect_equal(attr(combined.cat5, "abundance"), "mixed")
  out5 <- 
    PlotCatalogToPdf(combined.cat5, file = file.path(tempdir(), "test5.pdf"))
  expect_equal(out5, TRUE)
  expect_error(TransformCatalog(combined.cat5, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "abundance" attribute
  attr(cat10, "abundance") <- NULL
  combined.cat6 <- cbind(cat9, cat10)
  expect_equal(attr(combined.cat6, "abundance"), NULL)
  out6 <- 
    PlotCatalogToPdf(combined.cat6, file = file.path(tempdir(), "test6.pdf"))
  expect_equal(out6, TRUE)
  expect_error(TransformCatalog(combined.cat6, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # Test when performing cbind operation to catalogs with different "region"
  # attribute
  cat11 <- cat12 <- cat1
  
  # When catalogs have non-NULL different "region" attribute
  attr(cat11, "region") <- "exome"
  combined.cat7 <- cbind(cat1, cat11)
  expect_equal(attr(combined.cat7, "region"), "mixed")
  out7 <- 
    PlotCatalogToPdf(combined.cat7, file = file.path(tempdir(), "test7.pdf"))
  expect_equal(out7, TRUE)
  expect_error(TransformCatalog(combined.cat7, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When some catalog have NULL "region" attribute
  attr(cat11, "region") <- NULL
  combined.cat8 <- cbind(cat1, cat11)
  expect_equal(attr(combined.cat8, "region"), "mixed")
  out8 <- 
    PlotCatalogToPdf(combined.cat8, file = file.path(tempdir(), "test8.pdf"))
  expect_equal(out8, TRUE)
  expect_error(TransformCatalog(combined.cat8, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "region" attribute
  attr(cat12, "region") <- NULL
  combined.cat9 <- cbind(cat11, cat12)
  expect_equal(attr(combined.cat9, "region"), NULL)
  out9 <- 
    PlotCatalogToPdf(combined.cat9, file = file.path(tempdir(), "test9.pdf"))
  expect_equal(out9, TRUE)
  expect_error(TransformCatalog(combined.cat9, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  graphics.off()
  for (i in 0:9) {
    path <- file.path(tempdir(), paste0("test", i, ".pdf"))
    unlink(path)
  }
})

test_that("cbind method for DBS78Catalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  catalog.counts <- ReadCatalog("testdata/regress.cat.dbs.78.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)

  cat1 <- catalog.counts[, 1, drop = FALSE]
  cat2 <- catalog.counts[, 2, drop = FALSE]
  cat3 <- cbind(cat1, cat2)
  out <- PlotCatalogToPdf(cat3, file = file.path(tempdir(), "test0.pdf"))
  expect_equal(out, TRUE)
  
  # Test when performing cbind operation to catalogs with different "ref.genome"
  # attribute
  cat4 <- catalog.counts[, 3, drop = FALSE]
  
  # When some catalog has NULL "ref.genome" attribute
  attr(cat4, "ref.genome") <- NULL
  combined.cat <- cbind(cat1, cat4)
  expect_equal(attr(combined.cat, "ref.genome"), "mixed")
  out1 <- 
    PlotCatalogToPdf(combined.cat, file = file.path(tempdir(), "test1.pdf"))
  expect_equal(out1, TRUE)
  expect_error(TransformCatalog(combined.cat, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "ref.genome" attribute
  cat5 <- cat4
  combined.cat1 <- cbind(cat4, cat5)
  expect_equal(attr(combined.cat1, "ref.genome"), NULL)
  out2 <- 
    PlotCatalogToPdf(combined.cat1, file = file.path(tempdir(), "test2.pdf"))
  expect_equal(out2, TRUE)
  t.cat1 <- TransformCatalog(combined.cat1, target.ref.genome = "hg19",
                             target.catalog.type = "counts.signature")
  
  # When catalogs have non-NULL different "ref.genome" attribute
  cat6 <- cat1
  if (!"" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38")) {
    attr(cat6, "ref.genome") <-
      BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    combined.cat2 <- cbind(cat2, cat6)
    expect_equal(attr(combined.cat2, "ref.genome"), "mixed")
    out3 <- 
      PlotCatalogToPdf(combined.cat2, file = file.path(tempdir(), "test3.pdf"))
    expect_equal(out3, TRUE)
    unlink(file.path(tempdir(), "test3.pdf"))
    expect_error(TransformCatalog(combined.cat2, target.ref.genome = "hg19",
                                  target.catalog.type = "counts.signature"))
  }
  
  # Test when performing cbind operation to catalogs with different
  # "catalog.type" attribute
  cat7 <- cat8 <- TransformCatalog(cat1, target.catalog.type = "density")
  
  # When catalogs have non-NULL different "catalog.type" attribute
  expect_error(cbind(cat2, cat7))
  
  # When some catalog has NULL "catalog.type" attribute
  attr(cat7, "catalog.type") <- NULL
  expect_error(cbind(cat2, cat7))
  
  # When all catalogs have NULL "catalog.type" attribute
  attr(cat8, "catalog.type") <- NULL
  expect_error(cbind(cat7, cat8))
  expect_error(TransformCatalog(cat8, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # Test when performing cbind operation to catalogs with different "abundance"
  # attribute
  cat9 <- cat10 <- cat1
  
  # When catalogs have non-NULL different "abundance" attribute
  attr(cat9, "abundance") <- abundance.3bp.genome.unstranded.GRCh38
  combined.cat4 <- cbind(cat1, cat9)
  expect_equal(attr(combined.cat4, "abundance"), "mixed")
  out4 <- 
    PlotCatalogToPdf(combined.cat4, file = file.path(tempdir(), "test4.pdf"))
  expect_equal(out4, TRUE)
  expect_error(TransformCatalog(combined.cat4, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When some catalog have NULL "abundance" attribute
  attr(cat9, "abundance") <- NULL
  combined.cat5 <- cbind(cat1, cat9)
  expect_equal(attr(combined.cat5, "abundance"), "mixed")
  out5 <- 
    PlotCatalogToPdf(combined.cat5, file = file.path(tempdir(), "test5.pdf"))
  expect_equal(out5, TRUE)
  expect_error(TransformCatalog(combined.cat5, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "abundance" attribute
  attr(cat10, "abundance") <- NULL
  combined.cat6 <- cbind(cat9, cat10)
  expect_equal(attr(combined.cat6, "abundance"), NULL)
  out6 <- 
    PlotCatalogToPdf(combined.cat6, file = file.path(tempdir(), "test6.pdf"))
  expect_equal(out6, TRUE)
  expect_error(TransformCatalog(combined.cat6, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # Test when performing cbind operation to catalogs with different "region"
  # attribute
  cat11 <- cat12 <- cat1
  
  # When catalogs have non-NULL different "region" attribute
  attr(cat11, "region") <- "exome"
  combined.cat7 <- cbind(cat1, cat11)
  expect_equal(attr(combined.cat7, "region"), "mixed")
  out7 <- 
    PlotCatalogToPdf(combined.cat7, file = file.path(tempdir(), "test7.pdf"))
  expect_equal(out7, TRUE)
  expect_error(TransformCatalog(combined.cat7, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When some catalog have NULL "region" attribute
  attr(cat11, "region") <- NULL
  combined.cat8 <- cbind(cat1, cat11)
  expect_equal(attr(combined.cat8, "region"), "mixed")
  out8 <- 
    PlotCatalogToPdf(combined.cat8, file = file.path(tempdir(), "test8.pdf"))
  expect_equal(out8, TRUE)
  expect_error(TransformCatalog(combined.cat8, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "region" attribute
  attr(cat12, "region") <- NULL
  combined.cat9 <- cbind(cat11, cat12)
  expect_equal(attr(combined.cat9, "region"), NULL)
  out9 <- 
    PlotCatalogToPdf(combined.cat9, file = file.path(tempdir(), "test9.pdf"))
  expect_equal(out9, TRUE)
  expect_error(TransformCatalog(combined.cat9, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  graphics.off()
  for (i in 0:9) {
    path <- file.path(tempdir(), paste0("test", i, ".pdf"))
    unlink(path)
  }
})

test_that("cbind method for DBS144Catalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  catalog.counts <- ReadCatalog("testdata/regress.cat.dbs.144.csv",
                                ref.genome = "GRCh37",
                                region = "transcript", 
                                catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)

  cat1 <- catalog.counts[, 1, drop = FALSE]
  cat2 <- catalog.counts[, 2, drop = FALSE]
  cat3 <- cbind(cat1, cat2)
  out <- PlotCatalogToPdf(cat3, file = file.path(tempdir(), "test0.pdf"))
  expect_equal(out, TRUE)
  
  # Test when performing cbind operation to catalogs with different "ref.genome"
  # attribute
  cat4 <- catalog.counts[, 3, drop = FALSE]
  
  # When some catalog has NULL "ref.genome" attribute
  attr(cat4, "ref.genome") <- NULL
  combined.cat <- cbind(cat1, cat4)
  expect_equal(attr(combined.cat, "ref.genome"), "mixed")
  out1 <- 
    PlotCatalogToPdf(combined.cat, file = file.path(tempdir(), "test1.pdf"))
  expect_equal(out1, TRUE)
  expect_error(TransformCatalog(combined.cat, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "ref.genome" attribute
  cat5 <- cat4
  combined.cat1 <- cbind(cat4, cat5)
  expect_equal(attr(combined.cat1, "ref.genome"), NULL)
  out2 <- 
    PlotCatalogToPdf(combined.cat1, file = file.path(tempdir(), "test2.pdf"))
  expect_equal(out2, TRUE)
  t.cat1 <- TransformCatalog(combined.cat1, target.ref.genome = "hg19",
                             target.catalog.type = "counts.signature")
  
  # When catalogs have non-NULL different "ref.genome" attribute
  cat6 <- cat1
  if (!"" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38")) {
    attr(cat6, "ref.genome") <-
      BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    combined.cat2 <- cbind(cat2, cat6)
    expect_equal(attr(combined.cat2, "ref.genome"), "mixed")
    out3 <- 
      PlotCatalogToPdf(combined.cat2, file = file.path(tempdir(), "test3.pdf"))
    expect_equal(out3, TRUE)
    unlink(file.path(tempdir(), "test3.pdf"))
    expect_error(TransformCatalog(combined.cat2, target.ref.genome = "hg19",
                                  target.catalog.type = "counts.signature"))
  }
  
  # Test when performing cbind operation to catalogs with different
  # "catalog.type" attribute
  cat7 <- cat8 <- TransformCatalog(cat1, target.catalog.type = "density")
  
  # When catalogs have non-NULL different "catalog.type" attribute
  expect_error(cbind(cat2, cat7))
  
  # When some catalog has NULL "catalog.type" attribute
  attr(cat7, "catalog.type") <- NULL
  expect_error(cbind(cat2, cat7))
  
  # When all catalogs have NULL "catalog.type" attribute
  attr(cat8, "catalog.type") <- NULL
  expect_error(cbind(cat7, cat8))
  expect_error(TransformCatalog(cat8, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # Test when performing cbind operation to catalogs with different "abundance"
  # attribute
  cat9 <- cat10 <- cat1
  
  # When catalogs have non-NULL different "abundance" attribute
  attr(cat9, "abundance") <- abundance.3bp.genome.unstranded.GRCh38
  combined.cat4 <- cbind(cat1, cat9)
  expect_equal(attr(combined.cat4, "abundance"), "mixed")
  out4 <- 
    PlotCatalogToPdf(combined.cat4, file = file.path(tempdir(), "test4.pdf"))
  expect_equal(out4, TRUE)
  expect_error(TransformCatalog(combined.cat4, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When some catalog have NULL "abundance" attribute
  attr(cat9, "abundance") <- NULL
  combined.cat5 <- cbind(cat1, cat9)
  expect_equal(attr(combined.cat5, "abundance"), "mixed")
  out5 <- 
    PlotCatalogToPdf(combined.cat5, file = file.path(tempdir(), "test5.pdf"))
  expect_equal(out5, TRUE)
  expect_error(TransformCatalog(combined.cat5, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "abundance" attribute
  attr(cat10, "abundance") <- NULL
  combined.cat6 <- cbind(cat9, cat10)
  expect_equal(attr(combined.cat6, "abundance"), NULL)
  out6 <- 
    PlotCatalogToPdf(combined.cat6, file = file.path(tempdir(), "test6.pdf"))
  expect_equal(out6, TRUE)
  expect_error(TransformCatalog(combined.cat6, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # Test when performing cbind operation to catalogs with different "region"
  # attribute
  cat11 <- cat12 <- cat1
  
  # When catalogs have non-NULL different "region" attribute
  attr(cat11, "region") <- "exome"
  combined.cat7 <- cbind(cat1, cat11)
  expect_equal(attr(combined.cat7, "region"), "mixed")
  out7 <- 
    PlotCatalogToPdf(combined.cat7, file = file.path(tempdir(), "test7.pdf"))
  expect_equal(out7, TRUE)
  expect_error(TransformCatalog(combined.cat7, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When some catalog have NULL "region" attribute
  attr(cat11, "region") <- NULL
  combined.cat8 <- cbind(cat1, cat11)
  expect_equal(attr(combined.cat8, "region"), "mixed")
  out8 <- 
    PlotCatalogToPdf(combined.cat8, file = file.path(tempdir(), "test8.pdf"))
  expect_equal(out8, TRUE)
  expect_error(TransformCatalog(combined.cat8, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "region" attribute
  attr(cat12, "region") <- NULL
  combined.cat9 <- cbind(cat11, cat12)
  expect_equal(attr(combined.cat9, "region"), NULL)
  out9 <- 
    PlotCatalogToPdf(combined.cat9, file = file.path(tempdir(), "test9.pdf"))
  expect_equal(out9, TRUE)
  expect_error(TransformCatalog(combined.cat9, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  graphics.off()
  for (i in 0:9) {
    path <- file.path(tempdir(), paste0("test", i, ".pdf"))
    unlink(path)
  }
})

test_that("cbind method for DBS136Catalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  catalog.counts <- ReadCatalog("testdata/regress.cat.dbs.136.csv",
                                ref.genome = "GRCh37",
                                region = "genome",
                                catalog.type = "counts")
  colnames(catalog.counts) <- paste0("HepG2_", 1 : 4)

  cat1 <- catalog.counts[, 1, drop = FALSE]
  cat2 <- catalog.counts[, 2, drop = FALSE]
  cat3 <- cbind(cat1, cat2)
  out <- PlotCatalogToPdf(cat3, file = file.path(tempdir(), "test0.pdf"))
  expect_equal(out, TRUE)
  
  # Test when performing cbind operation to catalogs with different "ref.genome"
  # attribute
  cat4 <- catalog.counts[, 3, drop = FALSE]
  
  # When some catalog has NULL "ref.genome" attribute
  attr(cat4, "ref.genome") <- NULL
  combined.cat <- cbind(cat1, cat4)
  expect_equal(attr(combined.cat, "ref.genome"), "mixed")
  out1 <- 
    PlotCatalogToPdf(combined.cat, file = file.path(tempdir(), "test1.pdf"))
  expect_equal(out1, TRUE)
  expect_error(TransformCatalog(combined.cat, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "ref.genome" attribute
  cat5 <- cat4
  combined.cat1 <- cbind(cat4, cat5)
  expect_equal(attr(combined.cat1, "ref.genome"), NULL)
  out2 <- 
    PlotCatalogToPdf(combined.cat1, file = file.path(tempdir(), "test2.pdf"))
  expect_equal(out2, TRUE)
  t.cat1 <- TransformCatalog(combined.cat1, target.ref.genome = "hg19",
                             target.catalog.type = "counts.signature")
  
  # When catalogs have non-NULL different "ref.genome" attribute
  cat6 <- cat1
  if (!"" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38")) {
    attr(cat6, "ref.genome") <-
      BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    combined.cat2 <- cbind(cat2, cat6)
    expect_equal(attr(combined.cat2, "ref.genome"), "mixed")
    out3 <- 
      PlotCatalogToPdf(combined.cat2, file = file.path(tempdir(), "test3.pdf"))
    expect_equal(out3, TRUE)
    unlink(file.path(tempdir(), "test3.pdf"))
    expect_error(TransformCatalog(combined.cat2, target.ref.genome = "hg19",
                                  target.catalog.type = "counts.signature"))
  }
  
  # Test when performing cbind operation to catalogs with different
  # "catalog.type" attribute
  cat7 <- cat8 <- TransformCatalog(cat1, target.catalog.type = "density")
  
  # When catalogs have non-NULL different "catalog.type" attribute
  expect_error(cbind(cat2, cat7))
  
  # When some catalog has NULL "catalog.type" attribute
  attr(cat7, "catalog.type") <- NULL
  expect_error(cbind(cat2, cat7))
  
  # When all catalogs have NULL "catalog.type" attribute
  attr(cat8, "catalog.type") <- NULL
  expect_error(cbind(cat7, cat8))
  expect_error(TransformCatalog(cat8, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # Test when performing cbind operation to catalogs with different "abundance"
  # attribute
  cat9 <- cat10 <- cat1
  
  # When catalogs have non-NULL different "abundance" attribute
  attr(cat9, "abundance") <- abundance.3bp.genome.unstranded.GRCh38
  combined.cat4 <- cbind(cat1, cat9)
  expect_equal(attr(combined.cat4, "abundance"), "mixed")
  out4 <- 
    PlotCatalogToPdf(combined.cat4, file = file.path(tempdir(), "test4.pdf"))
  expect_equal(out4, TRUE)
  expect_error(TransformCatalog(combined.cat4, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When some catalog have NULL "abundance" attribute
  attr(cat9, "abundance") <- NULL
  combined.cat5 <- cbind(cat1, cat9)
  expect_equal(attr(combined.cat5, "abundance"), "mixed")
  out5 <- 
    PlotCatalogToPdf(combined.cat5, file = file.path(tempdir(), "test5.pdf"))
  expect_equal(out5, TRUE)
  expect_error(TransformCatalog(combined.cat5, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "abundance" attribute
  attr(cat10, "abundance") <- NULL
  combined.cat6 <- cbind(cat9, cat10)
  expect_equal(attr(combined.cat6, "abundance"), NULL)
  out6 <- 
    PlotCatalogToPdf(combined.cat6, file = file.path(tempdir(), "test6.pdf"))
  expect_equal(out6, TRUE)
  expect_error(TransformCatalog(combined.cat6, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # Test when performing cbind operation to catalogs with different "region"
  # attribute
  cat11 <- cat12 <- cat1
  
  # When catalogs have non-NULL different "region" attribute
  attr(cat11, "region") <- "exome"
  combined.cat7 <- cbind(cat1, cat11)
  expect_equal(attr(combined.cat7, "region"), "mixed")
  out7 <- 
    PlotCatalogToPdf(combined.cat7, file = file.path(tempdir(), "test7.pdf"))
  expect_equal(out7, TRUE)
  expect_error(TransformCatalog(combined.cat7, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When some catalog have NULL "region" attribute
  attr(cat11, "region") <- NULL
  combined.cat8 <- cbind(cat1, cat11)
  expect_equal(attr(combined.cat8, "region"), "mixed")
  out8 <- 
    PlotCatalogToPdf(combined.cat8, file = file.path(tempdir(), "test8.pdf"))
  expect_equal(out8, TRUE)
  expect_error(TransformCatalog(combined.cat8, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "region" attribute
  attr(cat12, "region") <- NULL
  combined.cat9 <- cbind(cat11, cat12)
  expect_equal(attr(combined.cat9, "region"), NULL)
  out9 <- 
    PlotCatalogToPdf(combined.cat9, file = file.path(tempdir(), "test9.pdf"))
  expect_equal(out9, TRUE)
  expect_error(TransformCatalog(combined.cat9, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  graphics.off()
  for (i in 0:9) {
    path <- file.path(tempdir(), paste0("test", i, ".pdf"))
    unlink(path)
  }
})

test_that("cbind method for IndelCatalog", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  catalog.counts <- ReadCatalog("testdata/BTSG_WGS_PCAWG.indels.csv",
                                ref.genome = "GRCh37",
                                region = "genome", catalog.type = "counts")

  cat1 <- catalog.counts[, 1, drop = FALSE]
  cat2 <- catalog.counts[, 2, drop = FALSE]
  cat3 <- cbind(cat1, cat2)
  out <- PlotCatalogToPdf(cat3, file = file.path(tempdir(), "test0.pdf"))
  expect_equal(out, TRUE)
  
  # Test when performing cbind operation to catalogs with different "ref.genome"
  # attribute
  cat4 <- catalog.counts[, 3, drop = FALSE]
  
  # When some catalog has NULL "ref.genome" attribute
  attr(cat4, "ref.genome") <- NULL
  combined.cat <- cbind(cat1, cat4)
  expect_equal(attr(combined.cat, "ref.genome"), "mixed")
  out1 <- 
    PlotCatalogToPdf(combined.cat, file = file.path(tempdir(), "test1.pdf"))
  expect_equal(out1, TRUE)
  expect_error(TransformCatalog(combined.cat, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "ref.genome" attribute
  cat5 <- cat4
  combined.cat1 <- cbind(cat4, cat5)
  expect_equal(attr(combined.cat1, "ref.genome"), NULL)
  out2 <- 
    PlotCatalogToPdf(combined.cat1, file = file.path(tempdir(), "test2.pdf"))
  expect_equal(out2, TRUE)
  expect_error(TransformCatalog(combined.cat1, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When catalogs have non-NULL different "ref.genome" attribute
  cat6 <- cat1
  if (!"" == system.file(package = "BSgenome.Hsapiens.UCSC.hg38")) {
    attr(cat6, "ref.genome") <-
      BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    combined.cat2 <- cbind(cat2, cat6)
    expect_equal(attr(combined.cat2, "ref.genome"), "mixed")
    out3 <- 
      PlotCatalogToPdf(combined.cat2, file = file.path(tempdir(), "test3.pdf"))
    expect_equal(out3, TRUE)
    unlink(file.path(tempdir(), "test3.pdf"))
    expect_error(TransformCatalog(combined.cat2, target.ref.genome = "hg19",
                                  target.catalog.type = "counts.signature"))
  }
  
  # Test when performing cbind operation to catalogs with different
  # "catalog.type" attribute
  tmp <- apply(cat1, MARGIN = 2, function(x) x / sum(x))
  cat7 <- cat8 <- as.catalog(tmp, ref.genome = "hg19", 
                             catalog.type = "counts.signature")
  
  # When catalogs have non-NULL different "catalog.type" attribute
  expect_error(cbind(cat2, cat7))
  
  # When some catalog has NULL "catalog.type" attribute
  attr(cat7, "catalog.type") <- NULL
  expect_error(cbind(cat2, cat7))
  
  # When all catalogs have NULL "catalog.type" attribute
  attr(cat8, "catalog.type") <- NULL
  expect_error(cbind(cat7, cat8))
  expect_error(TransformCatalog(cat8, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # Test when all catalogs have NULL "abundance" attribute
  cat9 <- cat10 <- cat1
  combined.cat6 <- cbind(cat9, cat10)
  expect_equal(attr(combined.cat6, "abundance"), NULL)
  out6 <- 
    PlotCatalogToPdf(combined.cat6, file = file.path(tempdir(), "test6.pdf"))
  expect_equal(out6, TRUE)
  expect_error(TransformCatalog(combined.cat6, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # Test when performing cbind operation to catalogs with different "region"
  # attribute
  cat11 <- cat12 <- cat1
  
  # When catalogs have non-NULL different "region" attribute
  attr(cat11, "region") <- "exome"
  combined.cat7 <- cbind(cat1, cat11)
  expect_equal(attr(combined.cat7, "region"), "mixed")
  out7 <- 
    PlotCatalogToPdf(combined.cat7, file = file.path(tempdir(), "test7.pdf"))
  expect_equal(out7, TRUE)
  expect_error(TransformCatalog(combined.cat7, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When some catalog have NULL "region" attribute
  attr(cat11, "region") <- NULL
  combined.cat8 <- cbind(cat1, cat11)
  expect_equal(attr(combined.cat8, "region"), "mixed")
  out8 <- 
    PlotCatalogToPdf(combined.cat8, file = file.path(tempdir(), "test8.pdf"))
  expect_equal(out8, TRUE)
  expect_error(TransformCatalog(combined.cat8, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  # When all catalogs have NULL "region" attribute
  attr(cat12, "region") <- NULL
  combined.cat9 <- cbind(cat11, cat12)
  expect_equal(attr(combined.cat9, "region"), NULL)
  out9 <- 
    PlotCatalogToPdf(combined.cat9, file = file.path(tempdir(), "test9.pdf"))
  expect_equal(out9, TRUE)
  expect_error(TransformCatalog(combined.cat9, target.ref.genome = "hg19",
                                target.catalog.type = "counts.signature"))
  
  graphics.off()
  for (i in 0:9) {
    path <- file.path(tempdir(), paste0("test", i, ".pdf"))
    unlink(path)
  }
})

