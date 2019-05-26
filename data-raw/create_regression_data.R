# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

files <- c("data-raw/VCF/HepG2_Cis_1_SNVresult_rmDup.vcf",
           "data-raw/VCF/HepG2_Cis_2_SNVresult_rmDup.vcf",
           "data-raw/VCF/HepG2_Cis_3_SNVresult_rmDup.vcf",
           "data-raw/VCF/HepG2_Cis_4_SNVresult_rmDup.vcf")

catalog <-
  StrelkaSBSVCFFilesToCatalog(files, genome = "GRCh37",
                              trans.ranges = trans.ranges.GRCh37)

WriteCatSBS96(catalog$catSBS96,
              path = "tests/testthat/testdata/regress.cat.sbs.96.csv")
WriteCatSBS192(catalog$catSBS192,
               path = "tests/testthat/testdata/regress.cat.sbs.192.csv")
WriteCatSBS1536(catalog$catSBS1536,
                path = "tests/testthat/testdata/regress.cat.sbs.1536.csv")
WriteCatDBS78(catalog$catDBS78,
              path = "tests/testthat/testdata/regress.cat.dbs.78.csv")
WriteCatDBS144(catalog$catDBS144,
               path = "tests/testthat/testdata/regress.cat.dbs.144.csv")
WriteCatDBS136(catalog$catDBS136,
               path = "tests/testthat/testdata/regress.cat.dbs.136.csv")

