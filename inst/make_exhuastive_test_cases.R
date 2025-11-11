# Make a nearly exhuastive set of test cases for ICAMS indel classification

tmpmock = read.csv("inst/annotations.unique.rows.csv") |>
  dplyr::select(
    -trans.start.pos,
    -trans.end.pos,
    -trans.strand,
    -trans.Ensembl.gene.ID,
    -trans.gene.symbol,
    -POS2,
    -bothstrand
  )
mock_vcf = dplyr::select(
  tmpmock,
  CHROM,
  POS,
  REF,
  ALT,
  Koh89.annotate.class,
  Koh476.annotate.class
)
write.csv(mock_vcf, "inst/exhaustive_mock_vcf.csv")

# Some tests

devtools::load_all()
avcf1 = AnnotateIDVCF(
  mock_vcf[28:32, ],
  "hg19",
  flag.mismatches = 1,
  explain_indels = T
)

library(magrittr)
check1 = function(pos) {
  # rowser()
  dplyr::filter(mock_vcf, POS == pos) -> err_row
  AnnotateIDVCF(
    err_row,
    "hg19",
    flag.mismatches = 1,
    explain_indels = TRUE
  )$annotated.vcf |>
    dplyr::select(
      -trans.start.pos,
      -trans.end.pos,
      -trans.strand,
      -trans.Ensembl.gene.ID,
      -trans.gene.symbol,
      -POS2,
      -bothstrand,
      -count
    ) %>%
    dplyr::relocate(Koh89.annotate.class, .after = Koh_89) %>%
    dplyr::relocate(Koh476.annotate.class, .after = Koh_476) %>%
    dplyr::relocate(COSMIC_83, .before = prev_COSMIC_83) -> vv
  View(vv, title = paste0("Indel at ", pos))
  return(invisible(vv))
}


devtools::load_all()


to_save = rbind(
  check1("56119277"),
  check1("47927645"),
  check1(171877227),
  check1("184762332"),
  check1("46904975"),
  check1(35132063)
)
write.csv(to_save, "inst/examples_where_COSMIC_and_Koh_differ.csv")
