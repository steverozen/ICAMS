test_that("Select and process Mutect VCF variants only from specific chromosomes", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  file <- "testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf"
  list.of.vcfs1 <- ReadAndSplitVCFs(file, variant.caller = "mutect")
  expect_equal(nrow(list.of.vcfs1$SBS[[1]]), 1738)
  
  # Only select variants that are in chromosomes 1 
  list.of.vcfs2 <- ReadAndSplitVCFs(file, variant.caller = "mutect",
                                    chr.names.to.process = "1")
  expect_equal(nrow(list.of.vcfs2$SBS[[1]]), 999)
  
  catalogs1 <-
    VCFsToCatalogs(files = file, ref.genome = "hg19", variant.caller = "mutect",
                   region = "genome")
  expect_equal(colSums(catalogs1$catSBS96), 1738, check.attributes = FALSE)
  
  # Only select variants that are in chromosomes 1 
  catalogs2 <-
    VCFsToCatalogs(files = file, ref.genome = "hg19", variant.caller = "mutect",
                   region = "genome", chr.names.to.process = "1")
  expect_equal(colSums(catalogs2$catSBS96), 999, check.attributes = FALSE)
  
  catalogs3 <-
    VCFsToCatalogsAndPlotToPdf(files = file, output.dir = tempdir(),
                               ref.genome = "hg19", variant.caller = "mutect",
                               region = "genome")
  expect_equal(colSums(catalogs3$catSBS96), 1738, check.attributes = FALSE)
  
  # Only select variants that are in chromosomes 1 
  catalogs4 <-
    VCFsToCatalogsAndPlotToPdf(files = file, output.dir = tempdir(),
                               ref.genome = "hg19", variant.caller = "mutect",
                               region = "genome",
                               chr.names.to.process = "1")
  expect_equal(colSums(catalogs4$catSBS96), 999, check.attributes = FALSE)
  
  catalogs5 <-
    VCFsToZipFile(files = file, zipfile = file.path(tempdir(), "test.zip"),
                  ref.genome = "hg19", variant.caller = "mutect",
                  region = "genome")
  expect_equal(colSums(catalogs5$catSBS96), 1738, check.attributes = FALSE)
  
  # Only select variants that are in chromosomes 1 
  catalogs6 <-
    VCFsToZipFile(files = file, zipfile = file.path(tempdir(), "test.zip"),
                  ref.genome = "hg19", variant.caller = "mutect",
                  region = "genome", chr.names.to.process = c("1"))
  expect_equal(colSums(catalogs6$catSBS96), 999, check.attributes = FALSE)
  
  names <- grep("\\.pdf$", list.files(tempdir()), value = TRUE)
  unlink(file.path(tempdir(), names))
})

test_that("Select and process Strelka SBS VCF variants only from specific chromosomes", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  file <- "testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf"
  list.of.vcfs1 <- ReadAndSplitVCFs(file, variant.caller = "strelka")
  expect_equal(nrow(list.of.vcfs1$SBS[[1]]), 770)
  
  # Only select variants that are in chromosomes 2
  list.of.vcfs2 <- ReadAndSplitVCFs(file, variant.caller = "strelka",
                                    chr.names.to.process = "2")
  expect_equal(nrow(list.of.vcfs2$SBS[[1]]), 0)
  
  catalogs1 <-
    VCFsToCatalogs(files = file, ref.genome = "hg19", variant.caller = "strelka",
                   region = "genome")
  expect_equal(colSums(catalogs1$catSBS96), 770, check.attributes = FALSE)
  
  # Only select variants that are in chromosomes 2
  catalogs2 <-
    VCFsToCatalogs(files = file, ref.genome = "hg19", variant.caller = "strelka",
                   region = "genome", chr.names.to.process = "2")
  expect_equal(colSums(catalogs2$catSBS96), 0, check.attributes = FALSE)
  
  catalogs3 <-
    VCFsToCatalogsAndPlotToPdf(files = file, output.dir = tempdir(),
                               ref.genome = "hg19", variant.caller = "strelka",
                               region = "genome")
  expect_equal(colSums(catalogs3$catSBS96), 770, check.attributes = FALSE)
  
  # Only select variants that are in chromosomes 2
  catalogs4 <-
    VCFsToCatalogsAndPlotToPdf(files = file, output.dir = tempdir(),
                               ref.genome = "hg19", variant.caller = "strelka",
                               region = "genome",
                               chr.names.to.process = "2")
  expect_equal(colSums(catalogs4$catSBS96), 0, check.attributes = FALSE)
  
  catalogs5 <-
    VCFsToZipFile(files = file, zipfile = file.path(tempdir(), "test.zip"),
                  ref.genome = "hg19", variant.caller = "strelka",
                  region = "genome")
  expect_equal(colSums(catalogs5$catSBS96), 770, check.attributes = FALSE)
  
  # Only select variants that are in chromosomes 2
  catalogs6 <-
    VCFsToZipFile(files = file, zipfile = file.path(tempdir(), "test.zip"),
                  ref.genome = "hg19", variant.caller = "strelka",
                  region = "genome", chr.names.to.process = c("2"))
  expect_equal(colSums(catalogs6$catSBS96), 0, check.attributes = FALSE)
  
  names <- grep("\\.pdf$", list.files(tempdir()), value = TRUE)
  unlink(file.path(tempdir(), names))
})

test_that("Select and process Strelka ID VCF variants only from specific chromosomes", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  file <- "testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf"
  list.of.vcfs1 <- ReadAndSplitVCFs(file, variant.caller = "strelka")
  expect_equal(nrow(list.of.vcfs1$ID[[1]]), 408)
  
  # Only select variants that are in chromosomes 2
  list.of.vcfs2 <- ReadAndSplitVCFs(file, variant.caller = "strelka",
                                    chr.names.to.process = "2")
  expect_equal(nrow(list.of.vcfs2$ID[[1]]), 32)
  
  catalogs1 <-
    VCFsToCatalogs(files = file, ref.genome = "hg19", variant.caller = "strelka",
                   region = "genome")
  expect_equal(colSums(catalogs1$catID), 408, check.attributes = FALSE)
  
  # Only select variants that are in chromosomes 2
  catalogs2 <-
    VCFsToCatalogs(files = file, ref.genome = "hg19", variant.caller = "strelka",
                   region = "genome", chr.names.to.process = "2")
  expect_equal(colSums(catalogs2$catID), 32, check.attributes = FALSE)
  
  catalogs3 <-
    VCFsToCatalogsAndPlotToPdf(files = file, output.dir = tempdir(),
                               ref.genome = "hg19", variant.caller = "strelka",
                               region = "genome")
  expect_equal(colSums(catalogs3$catID), 408, check.attributes = FALSE)
  
  # Only select variants that are in chromosomes 2
  catalogs4 <-
    VCFsToCatalogsAndPlotToPdf(files = file, output.dir = tempdir(),
                               ref.genome = "hg19", variant.caller = "strelka",
                               region = "genome",
                               chr.names.to.process = "2")
  expect_equal(colSums(catalogs4$catID), 32, check.attributes = FALSE)
  
  catalogs5 <-
    VCFsToZipFile(files = file, zipfile = file.path(tempdir(), "test.zip"),
                  ref.genome = "hg19", variant.caller = "strelka",
                  region = "genome")
  expect_equal(colSums(catalogs5$catID), 408, check.attributes = FALSE)
  
  # Only select variants that are in chromosomes 2
  catalogs6 <-
    VCFsToZipFile(files = file, zipfile = file.path(tempdir(), "test.zip"),
                  ref.genome = "hg19", variant.caller = "strelka",
                  region = "genome", chr.names.to.process = c("2"))
  expect_equal(colSums(catalogs6$catID), 32, check.attributes = FALSE)
  
  names <- grep("\\.pdf$", list.files(tempdir()), value = TRUE)
  unlink(file.path(tempdir(), names))
})

