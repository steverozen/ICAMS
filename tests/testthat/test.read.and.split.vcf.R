context("Reading and splitting VCFs")

test_that(
  "ReadListOfMutectVCFs",
  { foo <- ReadListOfMutectVCFs("VeryShortMutect.vcf")
    expect_equal(dim(foo[[1]]), c(1851, 12))
  } )

test_that(
  "ReadListOfStrelkaSNSVCFs applied to Mutect VCF error",
  {
    expect_error(ReadListOfStrelkaSNSVCFs("VeryShortMutect.vcf"),
                 "does not appear to a Strelka VCF",
                 fixed = TRUE)
  }
)
