context("Test deprecated function")

test_that("VCFs to catalogs functions", {
  rlang::with_options(lifecycle_verbosity = "warning", {
    skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
    stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
    
    file1 <- "testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf"
    expect_warning(
      {catalogs1 <- MutectVCFFilesToCatalog(file1, ref.genome = "hg19",
                              region = "genome")})
    
    file2 <- "testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf"
    expect_warning(
      {catalogs1 <- StrelkaSBSVCFFilesToCatalog(file2, ref.genome = "hg19",
                                             region = "genome")},
      regexp = 'Please use `VCFsToCatalogs(variant.caller = "strelka")` instead',
      fixed = TRUE)
    
    file3 <- "testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf"
    expect_warning(
      {catalogs1 <- StrelkaIDVCFFilesToCatalog(file3, ref.genome = "hg19",
                                            region = "genome")},
      regexp = 'Please use `VCFsToCatalogs(variant.caller = "strelka")` instead',
      fixed = TRUE)
  })
})

test_that("VCFs to catalogs and plot to Pdf functions", {
  rlang::with_options(lifecycle_verbosity = "warning", {
    skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
    stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
    
    output.files <- file.path(tempdir(), c("mutect", "strelka.sbs", "strelka.id"))
    file1 <- "testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf"
    expect_warning(
      {catalogs1 <- MutectVCFFilesToCatalogAndPlotToPdf(file1, ref.genome = "hg19",
                                          region = "genome",
                                          output.file = output.files[1])},
      regexp = 'Please use `VCFsToCatalogsAndPlotToPdf(variant.caller = "mutect")` instead',
      fixed = TRUE)
    
    file2 <- "testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf"
    expect_warning(
      {catalogs1 <- StrelkaSBSVCFFilesToCatalogAndPlotToPdf(file2, ref.genome = "hg19",
                                              region = "genome",
                                              output.file = output.files[2])},
      regexp = 'Please use `VCFsToCatalogsAndPlotToPdf(variant.caller = "strelka")` instead',
      fixed = TRUE)
    
    file3 <- "testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf"
    expect_warning(
      {catalogs1 <- StrelkaIDVCFFilesToCatalogAndPlotToPdf(file3, ref.genome = "hg19",
                                             region = "genome",
                                             output.file = output.files[3])},
      regexp = 'Please use `VCFsToCatalogsAndPlotToPdf(variant.caller = "strelka")` instead',
      fixed = TRUE)
   output.files <- list.files(path = tempdir(), pattern = ".pdf$", full.names = TRUE)
   sapply(output.files, FUN = unlink)
  })
})

test_that("VCFs to zip file functions", {
  rlang::with_options(lifecycle_verbosity = "warning", {
    skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
    stopifnot(requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5"))
    
    zipfile.names <- file.path(tempdir(), paste0("test", 1:3, ".zip"))
    dir1 <- "testdata/Mutect-GRCh37/"
    expect_warning(
      {catalogs1 <- MutectVCFFilesToZipFile(dir1, ref.genome = "hg19", zipfile = zipfile.names[1],
                              region = "genome")},
      regexp = 'Please use `VCFsToZipFile(variant.caller = "mutect")` instead',
      fixed = TRUE)
    
    dir2 <- "testdata/Strelka-SBS-GRCh37/"
    expect_warning(
      {catalogs1 <- StrelkaSBSVCFFilesToZipFile(dir2, ref.genome = "hg19", zipfile = zipfile.names[2],
                                  region = "genome")},
      regexp = 'Please use `VCFsToZipFile(variant.caller = "strelka")` instead',
      fixed = TRUE)
    
    dir3 <- "testdata/Strelka-ID-GRCh37/"
    expect_warning(
      {catalogs1 <- StrelkaIDVCFFilesToZipFile(dir3, ref.genome = "hg19", zipfile = zipfile.names[3],
                                 region = "genome")},
      regexp = 'Please use `VCFsToZipFile(variant.caller = "strelka")` instead',
      fixed = TRUE)
    sapply(zipfile.names, FUN = unlink)
  })
})

test_that("Read and split VCFs functions", {
  rlang::with_options(lifecycle_verbosity = "warning", {
    file1 <- "testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf"
    expect_warning(
        {split.vcfs1 <- ReadAndSplitMutectVCFs(file1)},
      regexp = 'Please use `ReadAndSplitVCFs(variant.caller = "mutect")` instead',
      fixed = TRUE)
    
    file2 <- "testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf"
    expect_warning(
      {split.vcfs2 <- ReadAndSplitStrelkaSBSVCFs(file2)},
      regexp = 'Please use `ReadAndSplitVCFs(variant.caller = "strelka")` instead',
      fixed = TRUE)
    
    file3 <- "testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf"
    expect_warning(
      {split.vcfs3 <- ReadStrelkaIDVCFs(file3)},
      regexp = 'Please use `ReadAndSplitVCFs(variant.caller = "strelka")` instead',
      fixed = TRUE)
  })
})
