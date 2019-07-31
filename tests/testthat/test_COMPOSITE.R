context("COMPOSITE")

test_that("ac.catalog COMPOSITE", {
  out <- as.catalog(numeric(1697), infer.rownames = TRUE)
  expect_equal(attr(out, "region", exact = TRUE), "unknown")
  expect_equal(class(out), c("COMPOSITECatalog", "matrix"))
})

test_that("WriteCatalog ReadCatalog COMPOSITE", {
  ctg <- as.catalog(numeric(1697), infer.rownames = TRUE)
  out.file <- tempfile(pattern = "1697", fileext = ".csv")
  WriteCatalog(ctg, file = out.file)
  ctg2 <- ReadCatalog(out.file)
  expect_equal(ctg, ctg2)
  
})
