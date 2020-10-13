context("Tests for internal functions reading and splitting vcfs")

test_that("Test ReadMutectVCFs and SplitListOfMutectVCFs", {
  files <- list.files(path = "testdata/Mutect-GRCh37", full.names = TRUE)
  vcfs <- ReadMutectVCFs(files)
  expect_equal(dim(vcfs[[1]]), c(1851,13))
  
  split.vcfs <- SplitListOfMutectVCFs(vcfs)
  expect_null(split.vcfs$discarded.variants$`Mutect.GRCh37.s2`)
})

test_that("Test ReadStrelkaSBSVCFs and SplitListOfStrelkaSBSVCFs", {
  files <- list.files(path = "testdata/Strelka-SBS-GRCh37/", full.names = TRUE)
  vcfs <- ReadStrelkaSBSVCFs(files)
  expect_equal(dim(vcfs[[1]]), c(798, 21))
  
  split.vcfs <- SplitListOfStrelkaSBSVCFs(vcfs)
  expect_null(split.vcfs$discarded.variants$`Mutect.GRCh37.s2`)
})

test_that("Test ReadAndSplitVCFs for Mutect VCFs", {
  file1 <- "testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf"
  split.vcfs1.1 <- ReadAndSplitMutectVCFs(file1)
  split.vcfs1.2 <- ReadAndSplitVCFs(file1, variant.caller = "mutect")
  expect_equal(split.vcfs1.1, split.vcfs1.2)
  
  file2 <- "testdata/Mutect-GRCh37/Mutect.GRCh37.s3.vcf"
  split.vcfs2.1 <- ReadAndSplitMutectVCFs(file2)
  split.vcfs2.2 <- ReadAndSplitVCFs(file2, variant.caller = "mutect")
  expect_equal(split.vcfs2.1, split.vcfs2.2)
  
  file3 <- "testdata/Mutect-GRCh37/Mutect.GRCh37.s4.vcf"
  split.vcfs3.1 <- ReadAndSplitMutectVCFs(file3)
  split.vcfs3.2 <- ReadAndSplitVCFs(file3, variant.caller = "mutect")
  expect_equal(split.vcfs3.1, split.vcfs3.2)
  
  file4 <- "testdata/Mutect-GRCh37/Mutect.GRCh37.s5.vcf"
  split.vcfs4.1 <- ReadAndSplitMutectVCFs(file4)
  split.vcfs4.2 <- ReadAndSplitVCFs(file4, variant.caller = "mutect")
  expect_equal(split.vcfs4.1, split.vcfs4.2)
  
  file5 <- list.files(path = "testdata/Mutect-GRCh37/", full.names = TRUE)
  split.vcfs5.1 <- ReadAndSplitMutectVCFs(file5)
  split.vcfs5.2 <- ReadAndSplitVCFs(file5, variant.caller = "mutect")
  expect_equal(split.vcfs5.1, split.vcfs5.2)
})

test_that("Test ReadAndSplitVCFs for Strelka SBS VCFs", {
  file1 <- "testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf"
  split.vcfs1.1 <- ReadAndSplitStrelkaSBSVCFs(file1)
  split.vcfs1.2 <- ReadAndSplitVCFs(file1, variant.caller = "strelka")
  split.vcfs1.2$ID <- NULL
  expect_equivalent(split.vcfs1.1, split.vcfs1.2)
  
  file2 <- "testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s3.vcf"
  split.vcfs2.1 <- ReadAndSplitStrelkaSBSVCFs(file2)
  split.vcfs2.2 <- ReadAndSplitVCFs(file2, variant.caller = "strelka")
  split.vcfs2.2$ID <- NULL
  expect_equivalent(split.vcfs2.1, split.vcfs2.2)
  
  file3 <- "testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s4.vcf"
  split.vcfs3.1 <- ReadAndSplitStrelkaSBSVCFs(file3)
  split.vcfs3.2 <- ReadAndSplitVCFs(file3, variant.caller = "strelka")
  split.vcfs3.2$ID <- NULL
  expect_equivalent(split.vcfs3.1, split.vcfs3.2)
  
  file4 <- "testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s5.vcf"
  split.vcfs4.1 <- ReadAndSplitStrelkaSBSVCFs(file4)
  split.vcfs4.2 <- ReadAndSplitVCFs(file4, variant.caller = "strelka")
  split.vcfs4.2$ID <- NULL
  expect_equivalent(split.vcfs4.1, split.vcfs4.2)
  
  file5 <- list.files(path = "testdata/Strelka-SBS-GRCh37/", full.names = TRUE)
  split.vcfs5.1 <- ReadAndSplitStrelkaSBSVCFs(file5)
  split.vcfs5.2 <- ReadAndSplitVCFs(file5, variant.caller = "strelka")
  split.vcfs5.2$ID <- NULL
  expect_equivalent(split.vcfs5.1, split.vcfs5.2)
})

test_that("Test ReadAndSplitVCFs for Strelka ID VCFs", {
  file1 <- "testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf"
  vcfs1.1 <- ReadStrelkaIDVCFs(file1)
  vcfs1.2 <- ReadAndSplitVCFs(file1, variant.caller = "strelka")
  expect_equal(vcfs1.1, vcfs1.2$ID)
  
  file2 <- "testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s3.vcf"
  vcfs2.1 <- ReadStrelkaIDVCFs(file2)
  vcfs2.2 <- ReadAndSplitVCFs(file2, variant.caller = "strelka")
  vcfs2.3 <- 
    suppressWarnings(lapply(vcfs2.1, 
                            FUN = CheckAndRemoveDiscardedVariants))
  vcfs2.4 <- lapply(vcfs2.3, FUN = "[[", 1)
  expect_equal(vcfs2.4, vcfs2.2$ID)
  
  file3 <- "testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s4.vcf"
  vcfs3.1 <- ReadStrelkaIDVCFs(file3)
  vcfs3.2 <- ReadAndSplitVCFs(file3, variant.caller = "strelka")
  expect_equal(vcfs3.1, vcfs3.2$ID)
  
  file4 <- "testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s5.vcf"
  vcfs4.1 <- ReadStrelkaIDVCFs(file4)
  vcfs4.2 <- ReadAndSplitVCFs(file4, variant.caller = "strelka")
  expect_equal(vcfs4.1, vcfs4.2$ID)
  
  file5 <- list.files(path = "testdata/Strelka-ID-GRCh37/", full.names = TRUE)
  vcfs5.1 <- ReadStrelkaIDVCFs(file5)
  vcfs5.2 <- ReadAndSplitVCFs(file5, variant.caller = "strelka")
  vcfs5.3 <- 
    suppressWarnings(lapply(vcfs5.1, 
                            FUN = CheckAndRemoveDiscardedVariants))
  vcfs5.4 <- lapply(vcfs5.3, FUN = "[[", 1)
  expect_equal(vcfs5.4, vcfs5.2$ID)
})