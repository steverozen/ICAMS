context("Reading and splitting VCFs")

test_that(
  "ReadStrelkaSBSVCFs",
  {
    vcf <- ReadStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh37.vcf")
    vcf1 <- ReadStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh38.vcf")
    expect_equal(dim(vcf[[1]]), c(798, 20))
    expect_equal(dim(vcf1[[1]]), c(1574, 12))
  } )

test_that(
  "ReadStrelkaIDVCFs",
  {
    vcf <- ReadStrelkaIDVCFs("testdata/Strelka.ID.GRCh37.vcf")
    vcf1 <- ReadStrelkaIDVCFs("testdata/Strelka.ID.GRCh38.vcf")
    expect_equal(dim(vcf[[1]]), c(408, 19))
    expect_equal(dim(vcf1[[1]]), c(1574, 11))
  } )

test_that(
  "ReadMutectVCFs",
  { foo <- ReadMutectVCFs("testdata/Mutect.GRCh37.vcf")
    expect_equal(dim(foo[[1]]), c(1851, 12))
  } )

test_that(
  "ReadStrelkaSBSVCFs applied to Mutect VCF error",
  {
    expect_error(ReadStrelkaSBSVCFs("testdata/Mutect.GRCh37.vcf"),
                 "does not appear to be a Strelka VCF",
                 fixed = TRUE)

    expect_error(ReadStrelkaSBSVCFs("testdata/Mutect.GRCh38.vcf"),
                 "does not appear to be a Strelka VCF",
                 fixed = TRUE)
  })

test_that(
  "ReadMutectVCFs applied to Strelka VCF error",
  {
    expect_error(ReadMutectVCFs("testdata/Strelka.SBS.GRCh37.vcf"),
                 "does not appear to be a Mutect VCF",
                 fixed = TRUE)

    expect_error(ReadMutectVCFs("testdata/Strelka.SBS.GRCh38.vcf"),
                 "does not appear to be a Mutect VCF",
                 fixed = TRUE)
  }
)
