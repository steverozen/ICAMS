context("Reading and splitting VCFs")

test_that(
  "ReadStrelkaSNSVCFs",
  {
    vcf <- ReadStrelkaSNSVCFs("testdata/Strelka.SNS.GRCh37.vcf")
    vcf1 <- ReadStrelkaSNSVCFs("testdata/Strelka.SNS.GRCh38.vcf")
    expect_equal(dim(vcf[[1]]), c(798, 20))
    expect_equal(dim(vcf1[[1]]), c(992, 12))
  } )

test_that(
  "ReadStrelkaIDVCFs",
  {
    vcf <- ReadStrelkaIDVCFs("testdata/Strelka.ID.GRCh37.vcf")
    vcf1 <- ReadStrelkaIDVCFs("testdata/Strelka.ID.GRCh38.vcf")
    expect_equal(dim(vcf[[1]]), c(408, 19))
    expect_equal(dim(vcf1[[1]]), c(406, 19))
  } )

test_that(
  "ReadMutectVCFs",
  { foo <- ReadMutectVCFs("testdata/Mutect.GRCh37.vcf")
    expect_equal(dim(foo[[1]]), c(1851, 12))
  } )

test_that(
  "ReadStrelkaSNSVCFs applied to Mutect VCF error",
  {
    expect_error(ReadStrelkaSNSVCFs("testdata/Mutect.GRCh37.vcf"),
                 "does not appear to be a Strelka VCF",
                 fixed = TRUE)
  })

test_that(
  "ReadMutectVCFs applied to Strelka VCF error",
  {
    expect_error(ReadMutectVCFs("testdata/Strelka.SNS.GRCh37.vcf"),
                 "does not appear to be a Mutect VCF",
                 fixed = TRUE)
  }
)
