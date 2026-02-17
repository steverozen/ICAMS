library(tidyr)
library(dplyr)

sample_id = "Colon::SP17905"

source("../Code_Liu_2025/code/read_annotated_vcf copy.R")

vv = read_annotated_vcf(sample_id)

vv %>%
  count(Koh_476, R, sort = TRUE, name = "newcount") -> vvv


cc = ICAMS::annot_vcf_to_476_catalog(vv)
