context("Reading and writing SigProfiler formatted Indel catalogs.")

test_that("Read SigProfiler", {
  expect_warning(
  sp.cat <-
    ReadCatalog(
      "testdata/sigProfiler_ID_signatures.csv",
      strict = FALSE)
  )
})

test_that("Write SigProfiler",{
  
  ICAMS.cat <- 
    ReadCatalog(
      "testdata/BTSG_WGS_PCAWG.indels.csv",
      strict = TRUE)
  
  ## Write Catalog in SigPro format
  WriteCatalogIndelSigPro(
    ICAMS.cat,
    "testdata/indels.sp.csv",
    strict = TRUE)
  
  ## Read Catalog in SigPro format, 
  ## obtain an ICAMS-formatted catalog matrix.
  expect_warning(
    ICAMS.cat.reread <-
      ReadCatalog(
        "testdata/indels.sp.csv",
        strict = FALSE)
  )
  
  ## Expect converted catalog to be equal to
  ## original catalog.
  expect_equal(ICAMS.cat,ICAMS.cat.reread)
  
  unlink("testdata/indels.sp.csv")
})
