# Make a nearly exhuastive set of test cases for ICAMS indel classification

testinput = "tests/testthat/testdata/long_test_vs_indelsiglib.csv"

if (FALSE) {
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
    seq.context,
    seq.context.width,
    Koh89.annotate.class,
    Koh476.annotate.class
  )
  readr::write_csv(mock_vcf, "inst/exhaustive_mock_vcf.csv")

  avcf1 = AnnotateIDVCF(
    mock_vcf[28:32, ],
    "hg19",
    flag.mismatches = 1,
    explain_indels = T
  )
}

# Some tests

check1 = function(pos) {
  mock_vcf = read.csv("inst/exhaustive_mock_vcf.csv")
  dplyr::filter(mock_vcf, POS == pos) -> err_row
  AnnotateIDVCF(
    err_row,
    "hg19",
    flag.mismatches = 1,
    explain_indels = TRUE
  )$annotated.vcf |>
    dplyr::relocate(Koh89.annotate.class, .after = Koh_89) |>
    dplyr::relocate(Koh476.annotate.class, .after = Koh_476) |>
    dplyr::relocate(COSMIC_83, .before = prev_COSMIC_83) -> vv
  View(vv, title = paste0("Indel at ", pos))
  return(invisible(vv))
}


test_indel_categorization = function() {
  mock_vcf = read.csv(testinput)
  retval1 = ICAMS:::categorize_indels_in_vcf(mock_vcf)
  # browser()
  cbind(mock_vcf, data.table::rbindlist(retval1, fill = TRUE)) -> bar

  bar |>
    dplyr::relocate(Koh89.annotate.class, .after = Koh_89) |>
    dplyr::relocate(Koh476.annotate.class, .after = Koh_476) |>
    dplyr::relocate(COSMIC_83, .before = prev_COSMIC_83) -> retval2
  # browser()
  diff1 = which(retval2$COSMIC_83 != retval2$prev_COSMIC_83)
  stopifnot(length(diff1) == 0)
  retval2$koh_orig_edited = make_koh_open_intervals(retval2)
  retval3 = dplyr::relocate(
    retval2,
    koh_orig_edited,
    .before = Koh89.annotate.class
  )
  retval3
}

make_koh_open_intervals = function(xx) {
  tmp = xx$Koh89.annotate.class
  tmp = gsub("Ins(2,):R(5,9)", "Ins(2,):R(5,)", x = tmp, fixed = TRUE) # Ins(C):R(7,)
  tmp = gsub("Ins(C):R(7,9)", "Ins(C):R(7,)", x = tmp, fixed = TRUE)
  tmp = gsub("[Del(T):R(8,9)]", "[Del(T):R(8,)]", x = tmp, fixed = TRUE)
  tmp = gsub("[Ins(T):R(8,9)]", "[Ins(T):R(8,)]", x = tmp, fixed = TRUE)
  tmp = gsub(
    "Del(3,):U(3,):R(3,)",
    "Del(3,):U(3,):R(3,)",
    x = tmp,
    fixed = TRUE
  )
  tmp = gsub(
    "Del(2,):U(1,2):R(5,9)",
    "Del(2,):U(1,2):R(5,)",
    x = tmp,
    fixed = TRUE
  )
  tmp = gsub(
    "Del(3,):U(3,):R(3,9)",
    "Del(3,):U(3,):R(3,)",
    x = tmp,
    fixed = TRUE
  )
  tmp = gsub(
    "Del(2,8):U(1,2):R(2,4)",
    "Del(2,):U(1,2):R(2,4)",
    x = tmp,
    fixed = TRUE
  )

  tmp
}

if (FALSE) {
  devtools::load_all()
  xx = test_indel_categorization()

  zz = which(xx$Koh_89 != xx$koh_orig_edited)
  length(zz)
  uu = (xx[zz, ])
  View(uu)
  # write.csv(xx[zz, ], file = "inst/our_koh89_ne_their_koh89.csv")
  # write.csv(xx, file = "inst/all_rows.csv")

  library(tidyr)

  xx %>% dplyr::filter(spacer_length == 0 & prime3_reps == 0 & L > 1) -> vv

  xx %>%
    dplyr::filter(stringr::str_detect(Koh89.annotate.class, ":M(\\d|\\()")) %>%
    dplyr::filter(spacer_length > 0) -> uu

  xx %>%
    dplyr::filter(!stringr::str_detect(Koh89.annotate.class, ":M(\\d|\\()")) %>%
    dplyr::filter(spacer_length > 0) -> uu # Only Ins...

  # examle single sample checking
  to_save = rbind(
    check1("56119277"),
    check1("47927645"),
    check1(171877227),
    check1("184762332"),
    check1("46904975"),
    check1(35132063)
  )
  write.csv(to_save, "inst/examples_where_COSMIC_and_Koh_differ.csv")
}
