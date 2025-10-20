setwd("C:/Users/steve/Documents/GitHub/ICAMS/tests/testthat")
source("../../inst/new_id_fns.R")
catalogs1 <- ICAMS::StrelkaIDVCFFilesToCatalog(
  list.files(path = "testdata/Strelka-ID-GRCh37/", full.names = TRUE),
  ref.genome = "hg19",
  region = "genome"
)
