test_that("AnnotateIDVCF_sample", {
  f1 = here::here("tests/testthat//testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf")
  vcf1 = ICAMS::ReadVCFs(f1, "mutect")[1]
  ivcf1 = dplyr::filter(vcf1[[1]], nchar(REF) != nchar(ALT))

  avcf1 = AnnotateIDVCF(
    ivcf1,
    "hg19",
    flag.mismatches = 0,
    explain_indels = 1
  )
  avcf1 = avcf1$annotated.vcf |>
    dplyr::select(U, R, Koh_89, Koh_476, COSMIC_83)
  expect_snapshot(avcf1)
  expect_snapshot(t(avcf1))
})
