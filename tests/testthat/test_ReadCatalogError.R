test_that("No file", {
  retval <- ReadCatalog("this.file.does.not.exist", stop.on.error = FALSE)
  expect_equal(
    retval,
    matrix(NA, nrow = 1, ncol = 1),
    check.attributes = FALSE
  )
})

test_that("Not a catalog file", {
  expect_warning({
    retval = ReadCatalog("test_ReadCatalogError.R", stop.on.error = FALSE)
  })
  expect_equal(
    retval,
    matrix(NA, nrow = 1, ncol = 1),
    check.attributes = FALSE
  )
})
