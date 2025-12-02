test_that("No file", {
  expect_error({retval <- ReadAndSplitVCFs("this.file.does.not.exist")})
})

test_that("Not a VCF file", {
  expect_error({retval <- ReadAndSplitVCFs("test_ReadVCFError.R",
                                           names.of.VCFs = "test_ReadVCFError.R")})
})