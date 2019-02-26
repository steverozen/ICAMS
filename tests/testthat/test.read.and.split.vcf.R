context("Reading and splitting VCFs")

test_that(
  "ReadMutectVCFs",
  { foo <- ReadMutectVCFs("VeryShortMutect.vcf")
    expect_equal(dim(foo[[1]]), c(1851, 12))
  } )

test_that(
  "ReadStrelkaSNSVCFs applied to Mutect VCF error",
  {
    expect_error(ReadStrelkaSNSVCFs("VeryShortMutect.vcf"),
                 "does not appear to a Strelka VCF",
                 fixed = TRUE)
  }
)
