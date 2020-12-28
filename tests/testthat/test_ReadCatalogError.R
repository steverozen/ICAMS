test_that("No file", {
  retval <- ReadCatalog("this.file.does.not.exist", stop.on.error = FALSE)
  expect_equal(
    retval,
    matrix(NA, nrow = 1, ncol = 1),
    check.attributes = FALSE)
})

test_that("Not a catalog file", {
  
  retval <- expect_warning(ReadCatalog("test_ReadCatalogError.R", stop.on.error = FALSE))
  expect_equal(
    retval,
    matrix(NA, nrow = 1, ncol = 1),
    check.attributes = FALSE)
})


if (FALSE) {
test_that("Bad row name", {
  retval <- ReadCatalog("testdata/regress.cat.sbs.96.bad.row.names.csv", stop.on.error = FALSE)
  expect_equal(
    retval,
    structure(NA, .Dim = c(1L, 1L), 
              error = 
                "File 'this.file.does.not.exist' does not exist or is non-readable. getwd()=='C:/Users/steve rozen/Documents/ICAMS'"))
})
}
