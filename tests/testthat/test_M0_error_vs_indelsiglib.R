test_that("test_M0_error_vs_indelsiglib", {
  cols_to_keep = c(
    "Sample",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "region",
    "mutID",
    "seq.context.width",
    "seq.context",
    "pos_shift",
    "ID.class",
    "dna.region",
    "ID166.class",
    "indelsigtool_Koh89"
  )

  mock_vcf = read.csv("testdata/M0_tests.tsv", sep = "\t") |>
    dplyr::select(all_of(cols_to_keep))
  retval1 = ICAMS:::categorize_indels_in_vcf(mock_vcf)
  cbind(mock_vcf, data.table::rbindlist(retval1, fill = TRUE)) -> xx
  zz = which(xx$Koh_89 != xx$indelsigtool_Koh89)
  # View(xx[zz, ] |> dplyr::select(Koh_89, indelsigtool_Koh89))
  expect_equal(length(zz), 0)
})
