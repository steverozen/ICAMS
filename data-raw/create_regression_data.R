# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

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

