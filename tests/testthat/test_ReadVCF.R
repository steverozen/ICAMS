context("Read VCFs")

test_that("Read Strelka mixed VCF", {
  file <- "testdata/Strelka.mixed.GRCh37.vcf"
  file1 <- "testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf"
  file2 <- "testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf"
  vcf <- ReadVCFs(files = file, variant.caller = "strelka")
  vcf1 <- ReadStrelkaSBSVCFs(files = file1)
  vcf2 <- ReadStrelkaIDVCFs(files = file2)
  rownames(vcf[[1]]) <- 1:nrow(vcf[[1]])
  expect_equal(vcf[[1]], dplyr::bind_rows(vcf1[[1]], vcf2[[1]]))
})

test_that(
  "Read Strelka SBS VCFs",
  {
    vcf <- ReadStrelkaSBSVCFs("testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf")
    vcf1 <- ReadStrelkaSBSVCFs("testdata/Strelka.SBS.GRCh38.vcf")
    vcf2 <- ReadStrelkaSBSVCFs("testdata/Strelka.SBS.GRCm38.vcf")
    expect_equal(dim(vcf[[1]]), c(798, 21))
    expect_equal(dim(vcf1[[1]]), c(1574, 13))
    expect_equal(dim(vcf2[[1]]), c(2870, 21))
    
    list <- ReadVCFs("testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf",
                     variant.caller = "strelka")
    list1 <- ReadVCFs("testdata/Strelka.SBS.GRCh38.vcf",
                      variant.caller = "strelka")
    list2 <- ReadVCFs("testdata/Strelka.SBS.GRCm38.vcf",
                      variant.caller = "strelka")
    expect_equal(vcf[[1]], list[[1]])
    expect_equal(vcf1[[1]], list1[[1]])
    expect_equal(vcf2[[1]], list2[[1]])
  } )

test_that(
  "Read Strelka ID VCFs",
  {
    vcf <- ReadStrelkaIDVCFs("testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf")
    vcf1 <- ReadStrelkaIDVCFs("testdata/Strelka.ID.GRCh38.vcf")
    vcf2 <- ReadStrelkaIDVCFs("testdata/Strelka.ID.GRCm38.vcf")
    expect_equal(dim(vcf[[1]]), c(408, 19))
    expect_equal(dim(vcf1[[1]]), c(1574, 11))
    expect_equal(dim(vcf2[[1]]), c(747, 19))
    
    list <- ReadVCFs("testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf",
                     variant.caller = "strelka")
    list1 <- ReadVCFs("testdata/Strelka.ID.GRCh38.vcf",
                     variant.caller = "strelka")
    list2 <- ReadVCFs("testdata/Strelka.ID.GRCm38.vcf",
                      variant.caller = "strelka")
    expect_equal(list[[1]], vcf[[1]])
    expect_equal(list1[[1]], vcf1[[1]])
    expect_equal(list2[[1]], vcf2[[1]])
  } )

test_that(
  "Read Mutect VCFs",
  { 
    vcf <- ReadMutectVCFs("testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf")
    vcf1 <- ReadMutectVCFs("testdata/Mutect.GRCh38.vcf")
    vcf2 <- expect_warning(ReadMutectVCFs("testdata/Mutect.GRCm38.vcf"))
    expect_equal(dim(vcf[[1]]), c(1851, 13))
    expect_equal(dim(vcf1[[1]]), c(1561, 13))
    expect_equal(dim(vcf2[[1]]), c(1895, 13))
    
    list <- ReadVCFs("testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf",
                     variant.caller = "mutect")
    list1 <- ReadVCFs("testdata/Mutect.GRCh38.vcf",
                      variant.caller = "mutect")
    list2 <- expect_warning(ReadVCFs("testdata/Mutect.GRCm38.vcf",
                                     variant.caller = "mutect",
                                     filter.status = NULL))
    expect_equal(list[[1]], vcf[[1]])
    expect_equal(list1[[1]], vcf1[[1]])
    expect_equal(list2[[1]], vcf2[[1]])
  } )

test_that(
  "ReadStrelkaSBSVCFs applied to Mutect VCF error",
  {
    expect_error(ReadStrelkaSBSVCFs("testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf"),
                 "does not appear to be a Strelka VCF",
                 fixed = TRUE)

    expect_error(ReadStrelkaSBSVCFs("testdata/Mutect.GRCh38.vcf"),
                 "does not appear to be a Strelka VCF",
                 fixed = TRUE)
    
    expect_error(ReadStrelkaSBSVCFs("testdata/Mutect.GRCm38.vcf"))
    
    expect_error(ReadVCFs("testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf",
                          variant.caller = "strelka"),
                 "does not appear to be a Strelka VCF",
                 fixed = TRUE)
    
    expect_error(ReadVCFs("testdata/Mutect.GRCh38.vcf",
                          variant.caller = "strelka"),
                 "does not appear to be a Strelka VCF",
                 fixed = TRUE)
    
    expect_error(ReadVCFs("testdata/Mutect.GRCm38.vcf",
                          variant.caller = "strelka"))
  })

test_that(
  "ReadStrelkaSBSVCFs applied to Strelka ID VCF error",
  {
    expect_error(ReadStrelkaSBSVCFs("testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf"),
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
    expect_error(ReadStrelkaIDVCFs("testdata/Mutect-GRCh37/Mutect.GRCh37.s1.vcf"),
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
    expect_error(ReadStrelkaIDVCFs("testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf"),
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
    expect_warning(ReadMutectVCFs("testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf"),
                   "does not appear to be a Mutect VCF",
                   fixed = TRUE)
    
    expect_warning(ReadMutectVCFs("testdata/Strelka.SBS.GRCh38.vcf"),
                   "does not appear to be a Mutect VCF",
                   fixed = TRUE)
    
    expect_warning(ReadMutectVCFs("testdata/Strelka.SBS.GRCm38.vcf"),
                   "does not appear to be a Mutect VCF",
                   fixed = TRUE)
    
    expect_warning(ReadMutectVCFs("testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf"),
                   "does not appear to be a Mutect VCF",
                   fixed = TRUE)
    
    expect_warning(ReadMutectVCFs("testdata/Strelka.ID.GRCh38.vcf"),
                 "does not appear to be a Mutect VCF",
                 fixed = TRUE)
    
    expect_warning(ReadMutectVCFs("testdata/Strelka.ID.GRCm38.vcf"),
                   "does not appear to be a Mutect VCF",
                   fixed = TRUE)
   
    expect_warning(ReadVCFs("testdata/Strelka-SBS-GRCh37/Strelka.SBS.GRCh37.s1.vcf",
                          variant.caller = "mutect"),
                 "does not appear to be a Mutect VCF",
                 fixed = TRUE)
    
    expect_warning(ReadVCFs("testdata/Strelka.SBS.GRCh38.vcf",
                          variant.caller = "mutect"),
                 "does not appear to be a Mutect VCF",
                 fixed = TRUE)
    
    expect_warning(ReadVCFs("testdata/Strelka.SBS.GRCm38.vcf",
                          variant.caller = "mutect"),
                 "does not appear to be a Mutect VCF",
                 fixed = TRUE)
    
    expect_warning(ReadVCFs("testdata/Strelka-ID-GRCh37/Strelka.ID.GRCh37.s1.vcf",
                          variant.caller = "mutect"),
                 "does not appear to be a Mutect VCF",
                 fixed = TRUE)
    
    expect_warning(ReadVCFs("testdata/Strelka.ID.GRCh38.vcf",
                          variant.caller = "mutect"),
                 "does not appear to be a Mutect VCF",
                 fixed = TRUE)
    
    expect_warning(ReadVCFs("testdata/Strelka.ID.GRCm38.vcf",
                            variant.caller = "mutect"),
                   "does not appear to be a Mutect VCF",
                   fixed = TRUE) 
  }
)
