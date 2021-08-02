context("Test deprecated function")

test_that("VCFs to catalogs functions", {
  rlang::with_options(lifecycle_verbosity = "warning", {
    file1 <- "testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf"
    catalogs1 <- expect_warning(
      MutectVCFFilesToCatalog(file1, ref.genome = "hg19", 
                              region = "genome"),
      regexp = 'Please use `VCFsToCatalogs(variant.caller = "mutect")` instead',
      fixed = TRUE)
    
    file2 <- "testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf"
    catalogs1 <- expect_warning(
      StrelkaSBSVCFFilesToCatalog(file2, ref.genome = "hg19", 
                                             region = "genome"),
      regexp = 'Please use `VCFsToCatalogs(variant.caller = "strelka")` instead',
      fixed = TRUE)
    
    file3 <- "testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf"
    catalogs1 <- expect_warning(
      StrelkaIDVCFFilesToCatalog(file3, ref.genome = "hg19", 
                                            region = "genome"),
      regexp = 'Please use `VCFsToCatalogs(variant.caller = "strelka")` instead',
      fixed = TRUE)
  })
})

test_that("VCFs to catalogs and plot to Pdf functions", {
  rlang::with_options(lifecycle_verbosity = "warning", {
    
    output.files <- file.path(tempdir(), c("mutect", "strelka.sbs", "strelka.id"))
    file1 <- "testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf"
    catalogs1 <- expect_warning(
      MutectVCFFilesToCatalogAndPlotToPdf(file1, ref.genome = "hg19", 
                                          region = "genome", 
                                          output.file = output.files[1]),
      regexp = 'Please use `VCFsToCatalogsAndPlotToPdf(variant.caller = "mutect")` instead',
      fixed = TRUE)
    
    file2 <- "testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf"
    catalogs1 <- expect_warning(
      StrelkaSBSVCFFilesToCatalogAndPlotToPdf(file2, ref.genome = "hg19", 
                                              region = "genome",
                                              output.file = output.files[2]),
      regexp = 'Please use `VCFsToCatalogsAndPlotToPdf(variant.caller = "strelka")` instead',
      fixed = TRUE)
    
    file3 <- "testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf"
    catalogs1 <- expect_warning(
      StrelkaIDVCFFilesToCatalogAndPlotToPdf(file3, ref.genome = "hg19", 
                                             region = "genome",
                                             output.file = output.files[3]),
      regexp = 'Please use `VCFsToCatalogsAndPlotToPdf(variant.caller = "strelka")` instead',
      fixed = TRUE)
   output.files <- list.files(path = tempdir(), pattern = ".pdf$", full.names = TRUE)
   sapply(output.files, FUN = unlink)
  })
})

test_that("VCFs to zip file functions", {
  rlang::with_options(lifecycle_verbosity = "warning", {
    zipfile.names <- file.path(tempdir(), paste0("test", 1:3, ".zip"))
    dir1 <- "testdata/Mutect-GRCh37/"
    catalogs1 <- expect_warning(
      MutectVCFFilesToZipFile(dir1, ref.genome = "hg19", zipfile = zipfile.names[1],
                              region = "genome"),
      regexp = 'Please use `VCFsToZipFile(variant.caller = "mutect")` instead',
      fixed = TRUE)
    
    dir2 <- "testdata/Strelka-SBS-GRCh37/"
    catalogs1 <- expect_warning(
      StrelkaSBSVCFFilesToZipFile(dir2, ref.genome = "hg19", zipfile = zipfile.names[2],
                                  region = "genome"),
      regexp = 'Please use `VCFsToZipFile(variant.caller = "strelka")` instead',
      fixed = TRUE)
    
    dir3 <- "testdata/Strelka-ID-GRCh37/"
    catalogs1 <- expect_warning(
      StrelkaIDVCFFilesToZipFile(dir3, ref.genome = "hg19", zipfile = zipfile.names[3],
                                 region = "genome"),
      regexp = 'Please use `VCFsToZipFile(variant.caller = "strelka")` instead',
      fixed = TRUE)
    sapply(zipfile.names, FUN = unlink)
  })
})

test_that("Read and split VCFs functions", {
  rlang::with_options(lifecycle_verbosity = "warning", {
    file1 <- "testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf"
    split.vcfs1 <- expect_warning(
        ReadAndSplitMutectVCFs(file1),
      regexp = 'Please use `ReadAndSplitVCFs(variant.caller = "mutect")` instead',
      fixed = TRUE)
    
    file2 <- "testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf"
    split.vcfs2 <- expect_warning(
      ReadAndSplitStrelkaSBSVCFs(file2),
      regexp = 'Please use `ReadAndSplitVCFs(variant.caller = "strelka")` instead',
      fixed = TRUE)
    
    file3 <- "testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf"
    split.vcfs3 <- expect_warning(
      ReadStrelkaIDVCFs(file3),
      regexp = 'Please use `ReadAndSplitVCFs(variant.caller = "strelka")` instead',
      fixed = TRUE)
  })
})
