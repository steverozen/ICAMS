context("Reading VCFs")

test_that(
  "ReadStrelkaSBSVCFs",
  {
    vcf <- ReadStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh37.vcf")
    vcf1 <- ReadStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh38.vcf")
    vcf2 <- ReadStrelkaSBSVCFs("testdata/Strelka.SBS.GRCm38.vcf")
    expect_equal(dim(vcf[[1]]), c(798, 20))
    expect_equal(dim(vcf1[[1]]), c(1574, 12))
    expect_equal(dim(vcf2[[1]]), c(2870, 20))
  } )

test_that(
  "ReadStrelkaIDVCFs",
  {
    vcf <- ReadStrelkaIDVCFs("testdata/Strelka.ID.GRCh37.vcf")
    vcf1 <- ReadStrelkaIDVCFs("testdata/Strelka.ID.GRCh38.vcf")
    vcf2 <- ReadStrelkaIDVCFs("testdata/Strelka.ID.GRCm38.vcf")
    expect_equal(dim(vcf[[1]]), c(408, 19))
    expect_equal(dim(vcf1[[1]]), c(1574, 11))
    expect_equal(dim(vcf2[[1]]), c(747, 19))
  } )

test_that(
  "ReadMutectVCFs",
  { vcf <- ReadMutectVCFs("testdata/Mutect-GRCh37/Mutect.GRCh37.vcf")
    vcf1 <- ReadMutectVCFs("testdata/Mutect.GRCh38.vcf")
    vcf2 <- expect_warning(ReadMutectVCFs("testdata/Mutect.GRCm38.vcf"))
    expect_equal(dim(vcf[[1]]), c(1851, 12))
    expect_equal(dim(vcf1[[1]]), c(1561, 12))
    expect_equal(dim(vcf2[[1]]), c(1895, 12))
  } )

test_that(
  "ReadStrelkaSBSVCFs applied to Mutect VCF error",
  {
    expect_error(ReadStrelkaSBSVCFs("testdata/Mutect-GRCh37/Mutect.GRCh37.vcf"),
                 "does not appear to be a Strelka VCF",
                 fixed = TRUE)

    expect_error(ReadStrelkaSBSVCFs("testdata/Mutect.GRCh38.vcf"),
                 "does not appear to be a Strelka VCF",
                 fixed = TRUE)
    
    expect_error(ReadStrelkaSBSVCFs("testdata/Mutect.GRCm38.vcf"))
  })

test_that(
  "ReadStrelkaSBSVCFs applied to Strelka ID VCF error",
  {
    expect_error(ReadStrelkaSBSVCFs("testdata/Strelka.ID.GRCh37.vcf"),
                 "does not appear to be a Strelka SBS VCF",
                 fixed = TRUE)
    
    expect_error(ReadStrelkaSBSVCFs("testdata/Strelka.ID.GRCh38.vcf"),
                 "does not appear to be a Strelka SBS VCF",
                 fixed = TRUE)
    
    expect_error(
      expect_warning(ReadStrelkaSBSVCFs("testdata/Strelka.ID.GRCm38.vcf")),
      "does not appear to be a Strelka SBS VCF",
      fixed = TRUE)
  })

test_that(
  "ReadStrelkaIDVCFs applied to Mutect VCF error",
  {
    expect_error(ReadStrelkaIDVCFs("testdata/Mutect-GRCh37/Mutect.GRCh37.vcf"),
                 "does not appear to be a Strelka VCF",
                 fixed = TRUE)
    
    expect_error(ReadStrelkaIDVCFs("testdata/Mutect.GRCh38.vcf"),
                 "does not appear to be a Strelka VCF",
                 fixed = TRUE)
    
    expect_error(ReadStrelkaIDVCFs("testdata/Mutect.GRCm38.vcf"))
  })

test_that(
  "ReadStrelkaIDVCFs applied to Strelka SBS VCF error",
  {
    expect_error(ReadStrelkaIDVCFs("testdata/Strelka.SBS.GRCh37.vcf"),
                 "does not appear to be a Strelka ID VCF",
                 fixed = TRUE)
    
    expect_error(ReadStrelkaIDVCFs("testdata/Strelka.SBS.GRCh38.vcf"),
                 "does not appear to be a Strelka ID VCF",
                 fixed = TRUE)
    
    expect_error(ReadStrelkaIDVCFs("testdata/Strelka.SBS.GRCm38.vcf"),
                 "does not appear to be a Strelka ID VCF",
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
    
    expect_error(ReadMutectVCFs("testdata/Strelka.SBS.GRCm38.vcf"),
                 "does not appear to be a Mutect VCF",
                 fixed = TRUE)
    
    expect_error(ReadMutectVCFs("testdata/Strelka.ID.GRCh37.vcf"),
                 "does not appear to be a Mutect VCF",
                 fixed = TRUE)
    
    expect_error(ReadMutectVCFs("testdata/Strelka.ID.GRCh38.vcf"),
                 "does not appear to be a Mutect VCF",
                 fixed = TRUE)
    
    expect_error(
      expect_warning(ReadMutectVCFs("testdata/Strelka.ID.GRCm38.vcf")),
      "does not appear to be a Mutect VCF",
      fixed = TRUE)
    
  }
)
