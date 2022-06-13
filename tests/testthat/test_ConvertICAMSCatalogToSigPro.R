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
   
   new.file <- file.path(tempdir(), "newfile.txt")
   ConvertCatalogToSigProfilerFormat(input.catalog = "testdata/regress.cat.sbs.96.csv",
                                     file = new.file)
   catalog <- ReadCatalog(file = new.file)
   expect_equal(input.catalog2, catalog)
   expect_true(tools::md5sum(file1) ==  tools::md5sum(new.file))
   
   unlink(c(file1, file2, new.file))
   
   # Testing converting catalog with only one column
   input.catalog3 <- "testdata/regress.cat.sbs.96.one.column.csv"
   file3 <- file.path(tempdir(), "test3.txt")
   ConvertICAMSCatalogToSigProSBS96(input.catalog = input.catalog3,
                                    file = file3)
   
   input.catalog4 <- ReadCatalog("testdata/regress.cat.sbs.96.one.column.csv")
   file4 <- file.path(tempdir(), "test4.txt")
   ConvertICAMSCatalogToSigProSBS96(input.catalog = input.catalog4,
                                    file = file4)
   expect_true(tools::md5sum(file3) ==  tools::md5sum(file4))
   
   file5 <- file.path(tempdir(), "test5.txt")
   input.catalog5 <- "testdata/regress.cat.sbs.96.one.column.csv"
   ConvertCatalogToSigProfilerFormat(input.catalog = input.catalog5,
                                     file = file5)
   catalog2 <- ReadCatalog(file = file5)
   expect_equal(catalog2, input.catalog4)
   expect_true(tools::md5sum(file3) ==  tools::md5sum(file5))
   
   unlink(c(file3, file4, file5))
})

test_that("Convert ICAMS Catalog to SigProfiler format SBS1536", {
  input.catalog1 <- "testdata/regress.cat.sbs.1536.csv"
  catalog1 <- ReadCatalog(input.catalog1)
  file1 <- file.path(tempdir(), "test1.txt")
  ConvertCatalogToSigProfilerFormat(input.catalog = "testdata/regress.cat.sbs.1536.csv",
                                    file = file1)
  catalog2 <- ReadCatalog(file = file1)
  expect_equal(catalog1, catalog2)
  
  unlink(file1)
  
  # Testing converting catalog with only one column
  input.catalog3 <- "testdata/regress.cat.sbs.1536.one.column.csv"
  catalog3 <- ReadCatalog(input.catalog3)
  file3 <- file.path(tempdir(), "test3.txt")
  ConvertCatalogToSigProfilerFormat(input.catalog = "testdata/regress.cat.sbs.1536.one.column.csv",
                                    file = file3)
  catalog4 <- ReadCatalog(file = file3)
  expect_equal(catalog3, catalog4)
  
  unlink(file3)
})

# TODO(Wuyang)
test_that("Convert ICAMS Catalog to SigProfiler format SBS192", {
  skip("Test for converting a SBS catalog is to be added.\n")
})


test_that("Convert ICAMS Catalog to SigProfiler format DBS78", {
   input.catalog1 <- "testdata/regress.cat.dbs.78.csv"
   catalog1 <- ReadCatalog(input.catalog1)
   file1 <- file.path(tempdir(), "test1.txt")
   ConvertCatalogToSigProfilerFormat(input.catalog = "testdata/regress.cat.dbs.78.csv",
                                     file = file1)
   catalog2 <- ReadCatalog(file = file1)
   expect_equal(catalog1, catalog2)
   
   unlink(file1)
   
   # Testing converting catalog with only one column
   input.catalog3 <- "testdata/regress.cat.dbs.78.one.column.csv"
   catalog3 <- ReadCatalog(input.catalog3)
   file3 <- file.path(tempdir(), "test3.txt")
   ConvertCatalogToSigProfilerFormat(input.catalog = "testdata/regress.cat.dbs.78.one.column.csv",
                                     file = file3)
   catalog4 <- ReadCatalog(file = file3)
   expect_equal(catalog3, catalog4)
   
   unlink(file3)
})

test_that("Convert ICAMS Catalog to SigProfiler format ID", {
   input.catalog1 <- "testdata/BTSG_WGS_PCAWG.indels.csv"
   catalog1 <- ReadCatalog(input.catalog1)
   file1 <- file.path(tempdir(), "test1.txt")
   ConvertCatalogToSigProfilerFormat(input.catalog = "testdata/BTSG_WGS_PCAWG.indels.csv",
                                     file = file1)
   catalog2 <- ReadCatalog(file = file1)
   expect_equal(catalog1, catalog2)
   
   unlink(file1)
   
   # Testing converting catalog with only one column
   input.catalog3 <- "testdata/BTSG_WGS_PCAWG.indels.one.column.csv"
   catalog3 <- ReadCatalog(input.catalog3)
   file3 <- file.path(tempdir(), "test3.txt")
   ConvertCatalogToSigProfilerFormat(input.catalog = "testdata/BTSG_WGS_PCAWG.indels.one.column.csv",
                                     file = file3)
   catalog4 <- ReadCatalog(file = file3)
   expect_equal(catalog3, catalog4)
   
   unlink(file3)
})

