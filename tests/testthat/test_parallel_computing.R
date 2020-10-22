context("Test parallel computing")

test_that("Test parallel computing for Strelka ID VCFs", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  files <- list.files(path = "testdata/Strelka-ID-GRCh37/", full.names = TRUE)

  if ("windows" == .Platform$OS.type) {
    num.of.cores <- 1
  } else {
    num.of.cores <- 2
  }

  list.of.ID.vcfs <-
    ReadVCFs(files, variant.caller = "strelka", num.of.cores = num.of.cores)
  catalogs1 <- VCFsToIDCatalogs(list.of.vcfs = list.of.ID.vcfs,
                                ref.genome = "hg19", region = "genome",
                                num.of.cores = num.of.cores)
  catalogs2 <- VCFsToCatalogs(files = files, ref.genome = "hg19",
                              variant.caller = "strelka",
                              num.of.cores = num.of.cores, region = "genome")
  expect_equal(catalogs1$catalog, catalogs2$catID)
  expect_equal(catalogs1$discarded.variants, catalogs2$discarded.variants)

  catalogs3 <- VCFsToIDCatalogs(list.of.vcfs = list.of.ID.vcfs,
                                ref.genome = "hg19", region = "genome",
                                num.of.cores = num.of.cores,
                                return.annotated.vcfs = TRUE)
  catalogs4 <- VCFsToCatalogs(files = files, ref.genome = "hg19",
                              variant.caller = "strelka",
                              num.of.cores = num.of.cores, region = "genome",
                              return.annotated.vcfs = TRUE)
  expect_equal(catalogs3$catalog, catalogs4$catID)
  expect_equal(catalogs3$discarded.variants, catalogs4$discarded.variants)
  expect_equal(catalogs3$annotated.vcfs, catalogs4$annotated.vcfs$ID)
})

test_that("Test parallel computing for Strelka SBS VCFs", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  files <- list.files(path = "testdata/Strelka-SBS-GRCh37/", full.names = TRUE)

  if ("windows" == .Platform$OS.type) {
    num.of.cores <- 1
  } else {
    num.of.cores <- 2
  }

  list.of.SBS.vcfs <-
    ReadVCFs(files, variant.caller = "strelka", num.of.cores = num.of.cores)

  split.vcfs <- SplitListOfVCFs(list.of.vcfs = list.of.SBS.vcfs,
                                num.of.cores = num.of.cores)

  catalogs1 <- VCFsToSBSCatalogs(list.of.SBS.vcfs = split.vcfs$SBS,
                                 ref.genome = "hg19",
                                 num.of.cores = num.of.cores,
                                 region = "genome")
  catalogs2 <- VCFsToCatalogs(files = files, ref.genome = "hg19",
                              variant.caller = "strelka",
                              num.of.cores = num.of.cores, region = "genome")
  expect_equal(catalogs1$catSBS96, catalogs2$catSBS96)

  catalogs3 <- VCFsToSBSCatalogs(list.of.SBS.vcfs = split.vcfs$SBS,
                                 ref.genome = "hg19", region = "genome",
                                 num.of.cores = num.of.cores,
                                 return.annotated.vcfs = TRUE)
  catalogs4 <- VCFsToCatalogs(files = files, ref.genome = "hg19",
                              variant.caller = "strelka",
                              num.of.cores = num.of.cores, region = "genome",
                              return.annotated.vcfs = TRUE)
  expect_equal(catalogs3$catSBS96, catalogs4$catSBS96)
  expect_equal(catalogs3$annotated.vcfs, catalogs4$annotated.vcfs$SBS)
})

test_that("Test parallel computing for Mutect VCFs", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  files <- list.files(path = "testdata/Mutect-GRCh37", full.names = TRUE)
  if ("windows" == .Platform$OS.type) {
    num.of.cores <- 1
  } else {
    num.of.cores <- 2
  }
  list.of.vcfs <-
    ReadVCFs(files, variant.caller = "mutect", num.of.cores = num.of.cores)
  split.vcfs <-
    SplitListOfVCFs(list.of.vcfs = list.of.vcfs, num.of.cores = num.of.cores)
  SBS.catalogs <- VCFsToSBSCatalogs(list.of.SBS.vcfs = split.vcfs$SBS,
                                    ref.genome = "hg19",
                                    num.of.cores = num.of.cores,
                                    region = "genome")

  DBS.catalogs <- VCFsToDBSCatalogs(list.of.DBS.vcfs = split.vcfs$DBS,
                                    ref.genome = "hg19",
                                    num.of.cores = num.of.cores,
                                    region = "genome")

  ID.catalogs <- VCFsToIDCatalogs(list.of.vcfs = split.vcfs$ID,
                                  ref.genome = "hg19",
                                  num.of.cores = num.of.cores,
                                  region = "genome")
  catalogs <- VCFsToCatalogs(files = files, ref.genome = "hg19",
                             variant.caller = "mutect",
                             num.of.cores = num.of.cores, region = "genome")
  expect_equal(SBS.catalogs$catSBS96, catalogs$catSBS96)
  expect_equal(DBS.catalogs$catDBS78, catalogs$catDBS78)
  expect_equal(ID.catalogs$catalog, catalogs$catID)
})
