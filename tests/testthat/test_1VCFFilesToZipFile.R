context("VCFFilestoZipfile functions")

test_that("MutectVCFFilesToZipFile function with no base.filename", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  dir <- "testdata/Mutect-GRCh37"
  out <- MutectVCFFilesToZipFile(dir, 
                                 zipfile = paste0(tempdir(), "/test.zip"), 
                                 ref.genome = "hg19", region = "genome")
  expect_type(out, "list")
  name <- grep("\\.zip$", list.files(tempdir()), value = TRUE)
  expect_equal(name, "test.zip")
  zip::unzip(zipfile = paste0(tempdir(), "/test.zip"), exdir = tempdir())
  name1 <- grep("\\.csv$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name1), 7)
  name2 <- grep("\\.pdf$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name2), 8)
  
  unlink(file.path(tempdir(), "test.zip"))
  unlink(file.path(tempdir(), name1))
  unlink(file.path(tempdir(), name2))
  graphics.off()
  unlink("testdata/Mutect-GRCh37/Rplots.pdf")
})

test_that("MutectVCFFilesToZipFile function with base.filename", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  dir <- "testdata/Mutect-GRCh37"
  out <- MutectVCFFilesToZipFile(dir, 
                                 zipfile = paste0(tempdir(), "/test.zip"), 
                                 ref.genome = "hg19",
                                 trans.ranges = trans.ranges.GRCh37,
                                 region = "genome",
                                 base.filename = "test")
  expect_type(out, "list")
  name <- grep("\\.zip$", list.files(tempdir()), value = TRUE)
  expect_equal(name, "test.zip")
  zip::unzip(zipfile = paste0(tempdir(), "/test.zip"), exdir = tempdir())
  name1 <- grep("\\.csv$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name1), 7)
  name2 <- grep("\\.pdf$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name2), 8)
  
  unlink(file.path(tempdir(), "test.zip"))
  unlink(file.path(tempdir(), name1))
  unlink(file.path(tempdir(), name2))
  graphics.off()
  unlink("testdata/Mutect-GRCh37/Rplots.pdf")
})

test_that("StrelkaSBSVCFFilesToZipFile function with no base.filename", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  dir <- "testdata/Strelka-SBS-GRCh37"
  out <- StrelkaSBSVCFFilesToZipFile(dir,
                                     zipfile = paste0(tempdir(), "/test.zip"), 
                                     ref.genome = "hg19",
                                     trans.ranges = trans.ranges.GRCh37,
                                     region = "genome")
  expect_type(out, "list")
  name <- grep("\\.zip$", list.files(tempdir()), value = TRUE)
  expect_equal(name, "test.zip")
  zip::unzip(zipfile = paste0(tempdir(), "/test.zip"), exdir = tempdir())
  name1 <- grep("\\.csv$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name1), 6)
  name2 <- grep("\\.pdf$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name2), 7)
  
  unlink(file.path(tempdir(), "test.zip"))
  unlink(file.path(tempdir(), name1))
  unlink(file.path(tempdir(), name2))
  graphics.off()
  unlink("testdata/Strelka-SBS-GRCh37/Rplots.pdf")
})

test_that("StrelkaSBSVCFFilesToZipFile function with base.filename", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  dir <- "testdata/Strelka-SBS-GRCh37"
  out <- StrelkaSBSVCFFilesToZipFile(dir,
                                     zipfile = paste0(tempdir(), "/test.zip"), 
                                     ref.genome = "hg19",
                                     trans.ranges = trans.ranges.GRCh37,
                                     region = "genome",
                                     base.filename = "test")
  expect_type(out, "list")
  name <- grep("\\.zip$", list.files(tempdir()), value = TRUE)
  expect_equal(name, "test.zip")
  zip::unzip(zipfile = paste0(tempdir(), "/test.zip"), exdir = tempdir())
  name1 <- grep("\\.csv$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name1), 6)
  name2 <- grep("\\.pdf$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name2), 7)
  
  unlink(file.path(tempdir(), "test.zip"))
  unlink(file.path(tempdir(), name1))
  unlink(file.path(tempdir(), name2))
  graphics.off()
  unlink("testdata/Strelka-SBS-GRCh37/Rplots.pdf")
})

test_that("StrelkaIDVCFFilesToZipFile function with no base.filename", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  dir <- "testdata/Strelka-ID-GRCh37"
  out <- StrelkaIDVCFFilesToZipFile(dir,
                                    zipfile = paste0(tempdir(), "/test.zip"),
                                    ref.genome = "hg19",
                                    region = "genome")
  expect_type(out, "list")
  name <- grep("\\.zip$", list.files(tempdir()), value = TRUE)
  expect_equal(name, "test.zip")
  
  zip::unzip(zipfile = paste0(tempdir(), "/test.zip"), exdir = tempdir())
  name1 <- grep("\\.csv$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name1), 1)
  name2 <- grep("\\.pdf$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name2), 1)
  
  unlink(file.path(tempdir(), "test.zip"))
  unlink(file.path(tempdir(), name1))
  unlink(file.path(tempdir(), name2))
  graphics.off()
  unlink("testdata/Strelka-ID-GRCh37/Rplots.pdf")
})

test_that("StrelkaIDVCFFilesToZipFile function with base.filename", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  dir <- "testdata/Strelka-ID-GRCh37"
  out <- StrelkaIDVCFFilesToZipFile(dir,
                                    zipfile = paste0(tempdir(), "/test.zip"),
                                    ref.genome = "hg19",
                                    region = "genome",
                                    base.filename = "test")
  expect_type(out, "list")
  name <- grep("\\.zip$", list.files(tempdir()), value = TRUE)
  expect_equal(name, "test.zip")
  
  zip::unzip(zipfile = paste0(tempdir(), "/test.zip"), exdir = tempdir())
  name1 <- grep("\\.csv$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name1), 1)
  name2 <- grep("\\.pdf$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name2), 1)
  
  unlink(file.path(tempdir(), "test.zip"))
  unlink(file.path(tempdir(), name1))
  unlink(file.path(tempdir(), name2))
  graphics.off()
  unlink("testdata/Strelka-ID-GRCh37/Rplots.pdf")
})

test_that("VCFsToZipFile function for Mutect VCFs with no base.filename", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  dir <- "testdata/Mutect-GRCh37"
  out <- VCFsToZipFile(dir, 
                       zipfile = paste0(tempdir(), "/test.zip"), 
                       ref.genome = "hg19", variant.caller = "mutect",
                       region = "genome")
  files <- list.files(path = dir, full.names = TRUE)
  out1 <- VCFsToZipFile(files = files, 
                        zipfile = paste0(tempdir(), "/test1.zip"), 
                        ref.genome = "hg19", variant.caller = "mutect",
                        region = "genome")
  expect_equal(out, out1)
  
  expect_type(out, "list")
  name <- grep("\\.zip$", list.files(tempdir()), value = TRUE)
  expect_equal(name, c("test.zip", "test1.zip"))
  zip::unzip(zipfile = paste0(tempdir(), "/test.zip"), exdir = tempdir())
  name1 <- grep("\\.csv$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name1), 7)
  name2 <- grep("\\.pdf$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name2), 8)
  
  unlink(file.path(tempdir(), "test.zip"))
  unlink(file.path(tempdir(), "test1.zip"))
  unlink(file.path(tempdir(), name1))
  unlink(file.path(tempdir(), name2))
  graphics.off()
  unlink("testdata/Mutect-GRCh37/Rplots.pdf")
})

test_that("VCFsToZipFile function for Mutect VCFs with base.filename", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  dir <- "testdata/Mutect-GRCh37"
  out <- VCFsToZipFile(dir, 
                       zipfile = paste0(tempdir(), "/test.zip"), 
                       ref.genome = "hg19",
                       region = "genome", variant.caller = "mutect",
                       base.filename = "test")
  files <- list.files(path = dir, full.names = TRUE)
  out1 <- VCFsToZipFile(files = files, 
                        zipfile = paste0(tempdir(), "/test1.zip"), 
                        ref.genome = "hg19",
                        region = "genome", variant.caller = "mutect",
                        base.filename = "test")
  expect_equal(out, out1)
  
  
  expect_type(out, "list")
  name <- grep("\\.zip$", list.files(tempdir()), value = TRUE)
  expect_equal(name, c("test.zip", "test1.zip"))
  zip::unzip(zipfile = paste0(tempdir(), "/test.zip"), exdir = tempdir())
  name1 <- grep("\\.csv$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name1), 7)
  name2 <- grep("\\.pdf$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name2), 8)
  
  unlink(file.path(tempdir(), "test.zip"))
  unlink(file.path(tempdir(), name1))
  unlink(file.path(tempdir(), name2))
  graphics.off()
  unlink("testdata/Mutect-GRCh37/Rplots.pdf")
})

test_that("VCFsToZipFile function for Strelka SBS VCFs with no base.filename", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  dir <- "testdata/Strelka-SBS-GRCh37"
  out <- VCFsToZipFile(dir,
                       zipfile = paste0(tempdir(), "/test.zip"), 
                       ref.genome = "hg19", variant.caller = "strelka",
                       region = "genome")
  
  files <- list.files(path = dir, full.names = TRUE)
  out1 <- VCFsToZipFile(files = files,
                        zipfile = paste0(tempdir(), "/test1.zip"), 
                        ref.genome = "hg19", variant.caller = "strelka",
                        region = "genome")
  expect_equal(out, out1)
  
  expect_type(out, "list")
  name <- grep("\\.zip$", list.files(tempdir()), value = TRUE)
  expect_equal(name, c("test.zip", "test1.zip"))
  zip::unzip(zipfile = paste0(tempdir(), "/test.zip"), exdir = tempdir())
  name1 <- grep("\\.csv$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name1), 7)
  name2 <- grep("\\.pdf$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name2), 7)
  
  unlink(file.path(tempdir(), "test.zip"))
  unlink(file.path(tempdir(), name1))
  unlink(file.path(tempdir(), name2))
  graphics.off()
  unlink("testdata/Strelka-SBS-GRCh37/Rplots.pdf")
})

test_that("VCFsToZipFile function for Strelka SBS VCFs with base.filename", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  dir <- "testdata/Strelka-SBS-GRCh37"
  out <- VCFsToZipFile(dir,
                       zipfile = paste0(tempdir(), "/test.zip"), 
                       ref.genome = "hg19",
                       region = "genome", variant.caller = "strelka",
                       base.filename = "test")
  
  files <- list.files(path = dir, full.names = TRUE)
  out1 <- VCFsToZipFile(files = files,
                        zipfile = paste0(tempdir(), "/test1.zip"), 
                        ref.genome = "hg19",
                        region = "genome", variant.caller = "strelka",
                        base.filename = "test")
  expect_equal(out, out1)
  
  expect_type(out, "list")
  name <- grep("\\.zip$", list.files(tempdir()), value = TRUE)
  expect_equal(name, c("test.zip", "test1.zip"))
  zip::unzip(zipfile = paste0(tempdir(), "/test.zip"), exdir = tempdir())
  name1 <- grep("\\.csv$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name1), 7)
  name2 <- grep("\\.pdf$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name2), 7)
  
  unlink(file.path(tempdir(), "test.zip"))
  unlink(file.path(tempdir(), name1))
  unlink(file.path(tempdir(), name2))
  graphics.off()
  unlink("testdata/Strelka-SBS-GRCh37/Rplots.pdf")
})

test_that("VCFsToZipFile function for Strelka ID VCFs with no base.filename", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  dir <- "testdata/Strelka-ID-GRCh37"
  out <- VCFsToZipFile(dir,
                       zipfile = paste0(tempdir(), "/test.zip"),
                       ref.genome = "hg19", variant.caller = "strelka",
                       region = "genome")
  
  files <- list.files(path = dir, full.names = TRUE)
  out1 <- VCFsToZipFile(files = files,
                        zipfile = paste0(tempdir(), "/test.zip"),
                        ref.genome = "hg19", variant.caller = "strelka",
                        region = "genome")
  expect_equal(out, out1)
  
  expect_type(out, "list")
  name <- grep("\\.zip$", list.files(tempdir()), value = TRUE)
  expect_equal(name, c("test.zip", "test1.zip"))
  
  zip::unzip(zipfile = paste0(tempdir(), "/test.zip"), exdir = tempdir())
  name1 <- grep("\\.csv$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name1), 7)
  name2 <- grep("\\.pdf$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name2), 1)
  
  unlink(file.path(tempdir(), "test.zip"))
  unlink(file.path(tempdir(), name1))
  unlink(file.path(tempdir(), name2))
  graphics.off()
  unlink("testdata/Strelka-ID-GRCh37/Rplots.pdf")
})

test_that("VCFsToZipFile function for Strelka ID VCFs with base.filename", {
  skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
  stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
  dir <- "testdata/Strelka-ID-GRCh37"
  out <- VCFsToZipFile(dir,
                       zipfile = paste0(tempdir(), "/test.zip"),
                       ref.genome = "hg19", variant.caller = "strelka",
                       region = "genome",
                       base.filename = "test")
  
  files <- list.files(path = dir, full.names = TRUE)
  out1 <- VCFsToZipFile(files = files,
                        zipfile = paste0(tempdir(), "/test1.zip"),
                        ref.genome = "hg19", variant.caller = "strelka",
                        region = "genome",
                        base.filename = "test")
  expect_equal(out, out1)
  
  expect_type(out, "list")
  name <- grep("\\.zip$", list.files(tempdir()), value = TRUE)
  expect_equal(name, c("test.zip", "test1.zip"))
  
  zip::unzip(zipfile = paste0(tempdir(), "/test.zip"), exdir = tempdir())
  name1 <- grep("\\.csv$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name1), 7)
  name2 <- grep("\\.pdf$", list.files(tempdir()), value = TRUE)
  expect_equal(length(name2), 1)
  
  unlink(file.path(tempdir(), "test.zip"))
  unlink(file.path(tempdir(), name1))
  unlink(file.path(tempdir(), name2))
  graphics.off()
  unlink("testdata/Strelka-ID-GRCh37/Rplots.pdf")
})