# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

# Create regression data in tests/testthat/testdata folder
files <- c("data-raw/VCF/HepG2_Cis_1_SNVresult_rmDup.vcf",
           "data-raw/VCF/HepG2_Cis_2_SNVresult_rmDup.vcf",
           "data-raw/VCF/HepG2_Cis_3_SNVresult_rmDup.vcf",
           "data-raw/VCF/HepG2_Cis_4_SNVresult_rmDup.vcf")

catalog <-
  StrelkaSBSVCFFilesToCatalog(files, ref.genome = "GRCh37",
                              trans.ranges = trans.ranges.GRCh37,
                              region = "genome")

WriteCatalog(catalog$catSBS96,
             file = "tests/testthat/testdata/regress.cat.sbs.96.csv")
WriteCatalog(catalog$catSBS192,
             file = "tests/testthat/testdata/regress.cat.sbs.192.csv")
WriteCatalog(catalog$catSBS1536,
             file = "tests/testthat/testdata/regress.cat.sbs.1536.csv")
WriteCatalog(catalog$catDBS78,
             file = "tests/testthat/testdata/regress.cat.dbs.78.csv")
WriteCatalog(catalog$catDBS144,
             file = "tests/testthat/testdata/regress.cat.dbs.144.csv")
WriteCatalog(catalog$catDBS136,
             file = "tests/testthat/testdata/regress.cat.dbs.136.csv")

# Create regression data in inst/extdata folder
file1 <- "data-raw/VCF/MCF10A_Carb_Low_cl2_Mutect.vcf"
catalog1 <- 
  MutectVCFFilesToCatalog(file1, ref.genome = "GRCh37",
                          trans.ranges = trans.ranges.GRCh37,
                          region = "genome")
WriteCatalog(catalog1$catSBS96,
             file = "inst/extdata/mutect.regress.cat.sbs.96.csv")
WriteCatalog(catalog1$catSBS192,
             file = "inst/extdata/mutect.regress.cat.sbs.192.csv")
WriteCatalog(catalog1$catSBS1536,
             file = "inst/extdata/mutect.regress.cat.sbs.1536.csv")
WriteCatalog(catalog1$catDBS78,
             file = "inst/extdata/mutect.regress.cat.dbs.78.csv")
WriteCatalog(catalog1$catDBS144,
             file = "inst/extdata/mutect.regress.cat.dbs.144.csv")
WriteCatalog(catalog1$catDBS136,
             file = "inst/extdata/mutect.regress.cat.dbs.136.csv")
WriteCatalog(catalog1$catID,
             file = "inst/extdata/mutect.regress.cat.indels.csv")

file2 <- "data-raw/VCF/MCF10A_Carb_Low_cl2_Strelka_SBS.vcf"
catalog2 <- 
  StrelkaSBSVCFFilesToCatalog(file2, ref.genome = "GRCh37",
                              trans.ranges = trans.ranges.GRCh37,
                              region = "genome")
WriteCatalog(catalog2$catSBS96,
             file = "inst/extdata/strelka.regress.cat.sbs.96.csv")
WriteCatalog(catalog2$catSBS192,
             file = "inst/extdata/strelka.regress.cat.sbs.192.csv")
WriteCatalog(catalog2$catSBS1536,
             file = "inst/extdata/strelka.regress.cat.sbs.1536.csv")
WriteCatalog(catalog2$catDBS78,
             file = "inst/extdata/strelka.regress.cat.dbs.78.csv")
WriteCatalog(catalog2$catDBS144,
             file = "inst/extdata/strelka.regress.cat.dbs.144.csv")
WriteCatalog(catalog2$catDBS136,
             file = "inst/extdata/strelka.regress.cat.dbs.136.csv")

file3 <- "data-raw/VCF/MCF10A_Carb_Low_cl2_Strelka_ID.vcf"
catalog3 <- StrelkaIDVCFFilesToCatalog(file3, ref.genome = "GRCh37", 
                                       region = "genome")
WriteCatalog(catalog3, file = "inst/extdata/strelka.regress.cat.indels.csv")
