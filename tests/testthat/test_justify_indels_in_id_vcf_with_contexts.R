test_that("justify_indels_in_id_vcf_with_contexts", {
  vcf_w_contexts = read.csv(
    "testdata/vcf_to_test_justify_indels_in_id_vcf_with_contexts.csv"
  )

  expect_snapshot(
    justify_indels_in_id_vcf_with_contexts(vcf_w_contexts, explain_indels = 2)
  )
})
