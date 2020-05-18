context("Convert ICAMS Catalog to SigProfiler format")

test_that("Convert ICAMS Catalog to SigProfiler format SBS96", {
   input.catalog1 <- "testdata/regress.cat.sbs.96.csv"
   file1 <- file.path(tempdir(), "test1.txt")
   ConvertICAMSCatalogToSigProSBS96(input.catalog = input.catalog1,
                                    file = file1)
   
   input.catalog2 <- ReadCatalog("testdata/regress.cat.sbs.96.csv")
   file2 <- file.path(tempdir(), "test2.txt")
   ConvertICAMSCatalogToSigProSBS96(input.catalog = input.catalog2,
                                    file = file2)
   expect_true(tools::md5sum(file1) ==  tools::md5sum(file2))
   
   unlink(c(file1, file2))
   
   # Testing converting catalog with only one column
   input.catalog3 <- "testdata/regress.cat.sbs.96.one.column.csv"
   file3 <- file.path(tempdir(), "test3.txt")
   ConvertICAMSCatalogToSigProSBS96(input.catalog = input.catalog3,
                                    file = file3)
   
   input.catalog4 <- "testdata/regress.cat.sbs.96.one.column.csv"
   file4 <- file.path(tempdir(), "test4.txt")
   ConvertICAMSCatalogToSigProSBS96(input.catalog = input.catalog4,
                                    file = file4)
   expect_true(tools::md5sum(file3) ==  tools::md5sum(file4))
   unlink(c(file3, file4))
})
