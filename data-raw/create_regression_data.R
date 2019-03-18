# Source this file from ICAMS top level directory.

cat(getwd(), "\n")

files <- c("data-raw/VCF/HepG2_Cis_1_SNVresult_rmDup.vcf",
           "data-raw/VCF/HepG2_Cis_2_SNVresult_rmDup.vcf",
           "data-raw/VCF/HepG2_Cis_3_SNVresult_rmDup.vcf",
           "data-raw/VCF/HepG2_Cis_4_SNVresult_rmDup.vcf")

catalog <-
  StrelkaSNSVCFFilesToCatalog(files, genome = "GRCh37",
                              trans.ranges = trans.ranges.GRCh37)

WriteCatSNS96(catalog$catSNS96,
              path = "tests/testthat/testdata/regress.cat.sns.96.csv")
WriteCatSNS192(catalog$catSNS192,
               path = "tests/testthat/testdata/regress.cat.sns.192.csv")
WriteCatSNS1536(catalog$catSNS1536,
                path = "tests/testthat/testdata/regress.cat.sns.1536.csv")
WriteCatDNS78(catalog$catDNS78,
              path = "tests/testthat/testdata/regress.cat.dns.78.csv")
WriteCatDNS144(catalog$catDNS144,
               path = "tests/testthat/testdata/regress.cat.dns.144.csv")
WriteCatDNS136(catalog$catDNS136,
               path = "tests/testthat/testdata/regress.cat.dns.136.csv")

