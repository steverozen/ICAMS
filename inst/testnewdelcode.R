library(ICAMS)
setwd("c:/Users/steve/Documents/GitHub/ICAMS/inst")


f1 = "c:/Users/steve/Documents/GitHub/ICAMS/tests/testthat/testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf"
vcf1 = ICAMS::ReadVCFs(f1, "mutect")[1]
ivcf1 = dplyr::filter(vcf1[[1]], nchar(REF) != nchar(ALT))

avcf1 = AnnotateIDVCF(
  ivcf1,
  "hg19",
  flag.mismatches = 0,
  explain_indels = F
)

avcf1 = avcf1$annotated.vcf
View(avcf1)

check1 = function(pos) {
  dplyr::filter(ivcf1, POS == pos) -> err_test
  AnnotateIDVCF(
    err_test,
    "hg19",
    flag.mismatches = 0,
    explain_indels = TRUE
  )$annotated.vcf |>
    View()
}
check1(146969286)
check1(102367067) # another Koh89 error


rdata = avcf1[, c(
  "CHROM",
  "POS",
  "REF",
  "ALT",
  "seq.context",
  "seq.context.width"
)]
ICAMS:::categorize_many_indels(rdata)

avcf1[, ICAMS:::gen_COSMIC_83_string(.SD), by = seq_len(nrow(avcf1))]

source("c:/Users/steve/Documents/GitHub/ICAMS/R/gen_koh_476_string.R")
avcf1[, gen_Koh_476_string(.SD), by = seq_len(nrow(avcf1))]

source("c:/Users/steve/Documents/GitHub/ICAMS/R/gen_koh_89_string.R")
avcf1[, xgen_Koh_89_string(.SD), by = seq_len(nrow(avcf1))]

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
