context("as.catalog")

test_that("as.catalog vector and no row names", {
 out <- as.catalog(numeric(1536), infer.rownames = TRUE)
 expect_equal(attr(out, "region", exact = TRUE), "unknown")
 expect_equal(class(out), c("SBS1536Catalog", "matrix"))
}
)

test_that("ac.catalog COMPOSITE", {
  out <- as.catalog(numeric(1697), infer.rownames = TRUE)
  expect_equal(attr(out, "region", exact = TRUE), "unknown")
  expect_equal(class(out), c("COMPOSITECatalog", "matrix"))
})


