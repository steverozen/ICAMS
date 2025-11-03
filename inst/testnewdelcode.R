library(ICAMS)
setwd("c:/Users/steve/Documents/GitHub/ICAMS/inst")


f1 = "c:/Users/steve/Documents/GitHub/ICAMS/tests/testthat/testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf"
vcf1 = ICAMS::ReadVCFs(f1, "mutect")[1]
ivcf1 = dplyr::filter(vcf1[[1]], nchar(REF) != nchar(ALT))

source("xcategorize_1_justified_indel.R")
avcf1 = AnnotateIDVCF(
  ivcf1,
  "hg19",
  flag.mismatches = 0,
  explain_indels = FALSE
)

avcf1 = avcf1$annotated.vcf
View(avcf1)


rdata = avcf1[, c(
  "CHROM",
  "POS",
  "REF",
  "ALT",
  "seq.context",
  "seq.context.width"
)]
ICAMS:::categorize_many_indels(rdata)

apply(avcf1, MARGIN = 1, FUN = ICAMS:::gen_COSMIC_83_string)

avcf1[, ICAMS:::gen_COSMIC_83_string(.SD), by = 1:nrow(avcf1)]

source("c:/Users/steve/Documents/GitHub/ICAMS/R/gen_koh_476_string.R")
avcf1[, xgen_Koh_476_string(.SD), by = 1:nrow(avcf1)]

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
