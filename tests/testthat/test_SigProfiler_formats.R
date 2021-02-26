context("Reading and writing SigProfiler formatted Indel catalogs.")

test_that("Read SigProfiler", {
  sp.cat <-
    ReadCatalog(
      "testdata/sigProfiler_ID_signatures.csv")
  expect_equal(dim(sp.cat), c(83, 17))
})

test_that("Write SigProfiler",{
  
  ICAMS.cat <- 
    ReadCatalog(
      "testdata/BTSG_WGS_PCAWG.indels.csv")
  
  ## Write Catalog in SigPro format
  WriteCatalogIndelSigPro(
    ICAMS.cat,
    "testdata/indels.sp.csv")
  
  ## Read Catalog in SigPro format, 
  ## obtain an ICAMS-formatted catalog matrix.
  ICAMS.cat.reread <-
    ReadCatalog(
      "testdata/indels.sp.csv")

  
  ## Expect converted catalog to be equal to
  ## original catalog.
  expect_equal(ICAMS.cat,ICAMS.cat.reread)
  
  unlink("testdata/indels.sp.csv")
})
