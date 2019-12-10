context("as.catalog")

test_that("as.catalog vector and no row names", {
 out <- as.catalog(numeric(1536), infer.rownames = TRUE)
 expect_equal(attr(out, "region", exact = TRUE), "unknown")
 expect_true(inherits(out, "SBS1536Catalog"))
}
)

test_that("as.catalog COMPOSITE", {
  out <- as.catalog(numeric(1697), infer.rownames = TRUE)
  expect_equal(attr(out, "region", exact = TRUE), "unknown")
  expect_true(inherits(out, "COMPOSITECatalog"))
})

test_that("as.catalog checks and reorders correct rownames", {
  mat <- matrix(data = 1:96, nrow = 96, ncol = 1,
                dimnames = list(catalog.row.order$SBS96))
  mat1 <- mat
  rownames(mat1)[1] <- "ACAU"
  expect_error(as.catalog(mat1))
  
  mat2 <- mat
  rownames(mat2)[2] <- "ACAA"
  expect_error(as.catalog(mat2))
  
  mat3 <- mat
  rownames(mat3)[1:2] <- c("ACAU", "ACCU")
  expect_error(as.catalog(mat3))
  
  mat4 <- mat[sort(rownames(mat)), , drop = FALSE]
  mat5 <- mat[sort(rownames(mat), decreasing = TRUE), , 
              drop = FALSE]
  cat1 <- as.catalog(mat4)
  cat2 <- as.catalog(mat5)
  expect_equal(cat1, cat2)
})


