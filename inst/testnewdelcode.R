library(ICAMS)
setwd("c:/Users/steve/Documents/GitHub/ICAMS/inst")


f1 = "c:/Users/steve/Documents/GitHub/ICAMS/tests/testthat/testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf"
vcf1 = ICAMS::ReadVCFs(f1, "mutect")[1]
ivcf1 = dplyr::filter(vcf1[[1]], nchar(REF) != nchar(ALT))

avcf1 = AnnotateIDVCF(
  ivcf1,
  "hg19",
  trans.ranges = NULL,
  flag.mismatches = 0,
  explain_indels = FALSE
)

avcf1 = avcf1$annotated.vcf

rdata = avcf1[, c(
  "CHROM",
  "POS",
  "REF",
  "ALT",
  "seq.context",
  "seq.context.width"
)]
ICAMS:::categorize_many_indels(rdata)


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
