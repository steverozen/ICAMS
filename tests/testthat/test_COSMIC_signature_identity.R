context("Testing whether signatures imported from COSMIC files are identical to signatures in PCAWG7 package.")

## This RDa file is saved from PCAWG7::signature object
## in PCAWG7 package version 0.0.3.9005.
load("testdata/PCAWGSig_0.0.3.9005.RDa")

test_that("Compare SBS1 in COSMIC-SBS96 format with SBS1 signature from PCAWG7 package in ICAMS-SBS96 format", {
  cosmic.SBS96 <-
    ReadCatalog(
      "testdata/COSMIC-Cat/COSMIC_SBS1_SBS96.csv",
	  catalog.type = "counts.signature")
  expect_equal(
    signif(cosmic.SBS96[,"SBS1_GRCh37"],3), 
	PCAWGSig$genome$SBS96[,"SBS1"])
})

test_that("Compare SBS1 in COSMIC-TSB192 format with SBS1 signature from PCAWG7 package in ICAMS-SBS192 format", {
  cosmic.TSB192 <-
    ReadCatalog(
      "testdata/COSMIC-Cat/COSMIC_SBS1_TSB192.csv",
	  catalog.type = "counts.signature")
  expect_equal(
    cosmic.TSB192[,"SBS1_GRCh37"], 
	PCAWGSig$genome$SBS192[,"SBS1"])
})

test_that("Compare DBS1 in COSMIC-DBS78 format with DBS1 signature from PCAWG7 package in ICAMS-DBS78 format", {
  cosmic.DBS78 <-
    ReadCatalog(
      "testdata/COSMIC-Cat/COSMIC_DBS1_DBS78.csv",
	  catalog.type = "counts.signature")
  expect_equal(
    signif(cosmic.DBS78[,"DBS1_GRCh37"],3),
	PCAWGSig$genome$DBS78[,"DBS1"])
})

test_that("Compare ID1 in COSMIC-ID83 format with ID1 signature from PCAWG7 package in ICAMS-ID83 format", {
  ## 13 channels not in ICAMS IndelCatalog format will be discarded.
  cosmic.ID83 <-
    ReadCatalog(
      "testdata/COSMIC-Cat/COSMIC_ID1_ID96.csv",
	  catalog.type = "counts.signature")
  expect_equal(
    cosmic.ID83[,"ID1_GRCh37"], 
	PCAWGSig$genome$ID[,"ID1"])
})
