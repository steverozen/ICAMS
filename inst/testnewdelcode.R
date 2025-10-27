library(ICAMS)
source("new_id_fns.R")
debug(xCanonicalize1ID)

files <- list.files(
  path = "c:/Users/steve/Documents/GitHub/ICAMS/tests/testthat/testdata/Mutect-GRCh37",
  full.names = TRUE
)
catalogs1 <- MutectVCFFilesToCatalog(
  files,
  ref.genome = "hg19",
  region = "genome"
)

# check function AnnotateIDVCF
# debug(AnnotateIDVCF)

f1 = "c:/Users/steve/Documents/GitHub/ICAMS/tests/testthat/testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf"
vcf1 = ICAMS::ReadVCFs(f1, "mutect")[1]
ivcf1 = dplyr::filter(vcf1[[1]], nchar(REF) != nchar(ALT))


source("annotate_ids_in_vcf.R")
debug(annotate_ids_in_vcf)
avcf1 = annotate_ids_in_vcf(
  ivcf1,
  "hg19",
  trans.ranges = NULL,
  flag.mismatches = 0
)


avcf1 = avcf1$annotated.vcf
