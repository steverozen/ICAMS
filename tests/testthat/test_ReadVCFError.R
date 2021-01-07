test_that("No file", {
  retval <- expect_error(ReadAndSplitVCFs("this.file.does.not.exist"))
})

test_that("Not a VCF file", {
  retval <- expect_error(ReadAndSplitVCFs("test_ReadVCFError.R", 
                                           names.of.VCFs = "test_ReadVCFError.R"))
})