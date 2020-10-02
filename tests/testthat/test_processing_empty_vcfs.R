context("Processing empty vcfs")

test_that("Processing empty Strelka SBS vcfs", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  
  file <- "testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s4.vcf"
  catalogs <- 
    StrelkaSBSVCFFilesToCatalog(file, ref.genome = "hg19", region = "genome")
  catalogs1 <- 
    StrelkaSBSVCFFilesToCatalog(file, ref.genome = "hg19", region = "genome",
                                return.annotated.vcfs = TRUE)
  expect_equivalent(colSums(catalogs$catSBS96), 0)
  
  file1 <- "testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s5.vcf"
  catalogs2 <- 
    StrelkaSBSVCFFilesToCatalog(file1, ref.genome = "hg19", region = "genome")
  catalogs3 <- 
    StrelkaSBSVCFFilesToCatalog(file1, ref.genome = "hg19", region = "genome",
                                return.annotated.vcfs = TRUE)
  expect_equivalent(colSums(catalogs2$catSBS96), 0)
  
  dir <- "testdata/Strelka-SBS-GRCh37/"
  catalogs4 <- 
    StrelkaSBSVCFFilesToZipFile(dir, 
                                zipfile = file.path(tempdir(), "strelka.sbs1.zip"),
                                ref.genome = "hg19", region = "genome")
  catalogs5 <- 
    StrelkaSBSVCFFilesToZipFile(dir, 
                                zipfile = file.path(tempdir(), "strelka.sbs2.zip"),
                                ref.genome = "hg19", region = "genome",
                                return.annotated.vcfs = TRUE)
  unlink(file.path(tempdir(), "strelka.sbs1.zip"))
  unlink(file.path(tempdir(), "strelka.sbs2.zip"))
})

test_that("Processing empty Strelka ID vcfs", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  
  file <- "testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s4.vcf"
  catalogs <- 
    StrelkaIDVCFFilesToCatalog(file, ref.genome = "hg19", region = "genome")
  catalogs1 <- 
    StrelkaIDVCFFilesToCatalog(file, ref.genome = "hg19", region = "genome",
                               return.annotated.vcfs = TRUE)
  expect_equivalent(colSums(catalogs$catalog), 0)
  
  file1 <- "testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s5.vcf"
  catalogs2 <- 
    StrelkaIDVCFFilesToCatalog(file1, ref.genome = "hg19", region = "genome")
  catalogs3 <- 
    StrelkaIDVCFFilesToCatalog(file1, ref.genome = "hg19", region = "genome",
                               return.annotated.vcfs = TRUE)
  expect_equivalent(colSums(catalogs2$catalog), 0)
  
  dir <- "testdata/Strelka-ID-GRCh37/"
  catalogs4 <- 
    StrelkaIDVCFFilesToZipFile(dir, 
                                zipfile = file.path(tempdir(), "strelka.id1.zip"),
                                ref.genome = "hg19", region = "genome")
  catalogs5 <- 
    StrelkaIDVCFFilesToZipFile(dir, 
                               zipfile = file.path(tempdir(), "strelka.id2.zip"),
                               ref.genome = "hg19", region = "genome",
                               return.annotated.vcfs = TRUE)
  unlink(file.path(tempdir(), "strelka.id1.zip"))
  unlink(file.path(tempdir(), "strelka.id2.zip"))
})

test_that("Processing empty Mutect vcfs", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  
  file <- "testdata/Mutect-GRCh37/Mutect.GRCh37.s4.vcf"
  catalogs <- 
    MutectVCFFilesToCatalog(file, ref.genome = "hg19", region = "genome")
  catalogs1 <- 
    MutectVCFFilesToCatalog(file, ref.genome = "hg19", region = "genome",
                            return.annotated.vcfs = TRUE)
  expect_equivalent(colSums(catalogs$catSBS96), 0)
  
  file1 <- "testdata/Mutect-GRCh37/Mutect.GRCh37.s5.vcf"
  catalogs2 <- 
    MutectVCFFilesToCatalog(file1, ref.genome = "hg19", region = "genome")
  catalogs3 <- 
    MutectVCFFilesToCatalog(file1, ref.genome = "hg19", region = "genome",
                            return.annotated.vcfs = TRUE)
  expect_equivalent(colSums(catalogs2$catSBS96), 0)
  
  dir <- "testdata/Mutect-GRCh37/"
  catalogs4 <- 
    MutectVCFFilesToZipFile(dir, 
                            zipfile = file.path(tempdir(), "mutect.1.zip"),
                            ref.genome = "hg19", region = "genome")
  catalogs5 <- 
    MutectVCFFilesToZipFile(dir, 
                            zipfile = file.path(tempdir(), "mutect.2.zip"),
                            ref.genome = "hg19", region = "genome",
                            return.annotated.vcfs = TRUE)
  unlink(file.path(tempdir(), "mutect.1.zip"))
  unlink(file.path(tempdir(), "mutect.2.zip"))
  graphics.off()
  unlink("Rplots.pdf")
})