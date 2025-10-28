library(ICAMS)
setwd("c:/Users/steve/Documents/GitHub/ICAMS/inst")


f1 = "c:/Users/steve/Documents/GitHub/ICAMS/tests/testthat/testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf"
vcf1 = ICAMS::ReadVCFs(f1, "mutect")[1]
ivcf1 = dplyr::filter(vcf1[[1]], nchar(REF) != nchar(ALT))

source("annotate_ids_in_vcf.R")
source("categorize_many_indels.R")
source("gen_COSMIC_83_string.R")
source("justify_indel.R")
source("xCanonicalize1ID.R")
source("categorize_1_justified_indel.R")
avcf1 = annotate_ids_in_vcf(
  ivcf1,
  "hg19",
  trans.ranges = NULL,
  flag.mismatches = 0,
  explain_indels = TRUE
)

avcf1 = avcf1$annotated.vcf

## Older, end-to-end tests

files <- list.files(
  path = "c:/Users/steve/Documents/GitHub/ICAMS/tests/testthat/testdata/Mutect-GRCh37",
  full.names = TRUE
)
catalogs1 <- MutectVCFFilesToCatalog(
  files,
  ref.genome = "hg19",
  region = "genome",
  return.annotated.vcfs = TRUE
)

files <- list.files(
  path = "c:/Users/steve/Documents/GitHub/ICAMS/tests/testthat/testdata/Strelka-ID-GRCh37/",
  full.names = TRUE
)
catalogs1 <- StrelkaIDVCFFilesToCatalog(
  files,
  ref.genome = "hg19",
  region = "genome"
)
