context("Test IsICAMSCatalog function")

test_that("IsICAMSCatalog function", {
dir <- c(system.file("extdata/Mutect-vcf",
                     package = "ICAMS"))
catalogs <-
  VCFsToZipFile(dir,
                zipfile = file.path(tempdir(), "test.zip"),
                ref.genome = "hg19",
                variant.caller = "mutect",
                region = "genome",
                base.filename = "Mutect")
cat.sbs96 <- catalogs$catSBS96
dim(cat.sbs96)
expect_true(IsICAMSCatalog(cat.sbs96))

# Test subsetting the catalog
test <- cat.sbs96[, 1, drop = FALSE]
dim(test)
expect_true(IsICAMSCatalog(test)) 

# Test a matrix with incorrect number of rows
mat1 <- matrix(data = 1, nrow = 97)
expect_false(IsICAMSCatalog(mat1)) 

# Test a matrix with correct number of rows but incorrect class
mat2 <- matrix(data = 1, nrow = 96)
expect_false(IsICAMSCatalog(mat2)) 

# Test a matrix with correct number of rows and class, but incorrect
# rownames
mat3 <- matrix(data = 1, nrow = 96)
mat3 <- as.catalog(mat3, infer.rownames = TRUE)
rownames(mat3) <- paste0("test", 1:96)

result <- 
  expect_message(IsICAMSCatalog(mat3), 
                 "The rownames of the input object do not match the catalog row order used in ICAMS exactly.") 
expect_false(result)
})