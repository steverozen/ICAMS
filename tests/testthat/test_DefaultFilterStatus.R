context("filter.status argument in reading VCFs")

test_that("Reading Strelka VCF", {
  file1 <- "testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s6.vcf"
  # User must specify the value of filter.status explicitly when variant.caller
  # is "unknown"
  vcf1 <- expect_error(ReadVCFs(file = file1))
  
  vcf2 <- ReadVCFs(file = file1, filter.status = "PASS")
  expect_equal(nrow(vcf2[[1]]), 3)
  
  # When variant.caller = "strelka", the filter.status will be "PASS" by default
  vcf3 <- ReadVCFs(file = file1, variant.caller = "strelka")
  expect_equal(nrow(vcf3[[1]]), 3)
  
  # When filter.status = NULL the code doesn't care if there is no FILTER column
  # All variants will be retained
  vcf4 <- ReadVCFs(file = file1, filter.status = NULL)
  expect_equal(nrow(vcf4[[1]]), 5)
  
  # When variant.caller is unrecognized, throw an error
  vcf5 <- expect_error(ReadVCFs(file = file1, variant.caller = "test"))
  
})

test_that("Reading Mutect VCF", {
  file1 <- "testdata/Mutect-GRCh37/Mutect.GRCh37.s3.vcf"
  # User must specify the value of filter.status explicitly when variant.caller
  # is "unknown"
  vcf1 <- expect_error(ReadVCFs(file = file1))
  
  vcf2 <- ReadVCFs(file = file1, filter.status = "PASS")
  expect_equal(nrow(vcf2[[1]]), 851)
  
  # When variant.caller = "strelka", the filter.status will be "PASS" by default
  vcf3 <- ReadVCFs(file = file1, variant.caller = "mutect")
  expect_equal(nrow(vcf3[[1]]), 851)
  
  # When filter.status = NULL the code doesn't care if there is no FILTER column
  # All variants will be retained
  vcf4 <- ReadVCFs(file = file1, filter.status = NULL)
  expect_equal(nrow(vcf4[[1]]), 851)
  
  # When variant.caller is unrecognized, throw an error
  vcf5 <- expect_error(ReadVCFs(file = file1, variant.caller = "test"))
})

test_that("VCF with no 'FILTER' column", {
  file1 <- "testdata/Strelka.SBS.GRCh38.no.FILTER.column.vcf"
  # User must specify the value of filter.status explicitly when variant.caller
  # is "unknown"
  vcf1 <- expect_error(ReadVCFs(file = file1))
  
  # When there is no FILTER column, generate a warning and all variants 
  # will be retained
  vcf2 <- expect_warning(ReadVCFs(file = file1, filter.status = "PASS"))
  expect_equal(nrow(vcf2[[1]]), 7)
  
  vcf3 <- expect_warning(ReadVCFs(file = file1, variant.caller = "strelka"))
  expect_equal(nrow(vcf3[[1]]), 7)
  
  # When filter.status = NULL the code doesn't care if there is no FILTER column
  vcf4 <- ReadVCFs(file = file1, filter.status = NULL)
  expect_equal(nrow(vcf4[[1]]), 7)
})