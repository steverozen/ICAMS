test_that("justify_indels_in_id_vcf_with_contexts", {
  vcf_no_contexts = read.csv(
    "testdata/vcf_to_test_justify_indels_in_id_vcf_with_contexts.csv"
  ) |>
    dplyr::select(-seq.context.width, -seq.context)

  expect_snapshot(
    justify_id_vcf(vcf_no_contexts, ref.genome = "hg19", explain_indels = 2)
  )
})
