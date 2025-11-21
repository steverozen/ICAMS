# tmptest_justify_id_vcf

if (FALSE) {
  # for setup
  xx = ICAMS::ReadVCFs(
    "tests/testthat/testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf",
    variant.caller = "strelka"
  )[[1]] |>
    dplyr::select(1:5)

  devtools::load_all()

  debug(ICAMS:::justify_indels_in_id_vcf_with_contexts)
  justify_id_vcf(xx, ref.genome = "hg19")
}

# Now we have a pseudo vcf with contexts that we can edit
vcf_w_contexts = read.csv(
  "inst/vcf_to_test_justify_indels_in_id_vcf_with_contexts.csv"
)

justify_indels_in_id_vcf_with_contexts(vcf_w_contexts, explain_indels = 2)
