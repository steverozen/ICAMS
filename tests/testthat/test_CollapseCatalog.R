context("CollapseCatalog")

test_that("Collapse1536CatalogTo96 function is working properly", {
  cat.SBS1536 <-
    ReadCatalog("testdata/regress.cat.sbs.1536.csv", ref.genome = "GRCh37",
                region = "genome", catalog.type = "counts")
  x1 <- Collapse1536CatalogTo96(cat.SBS1536)
  expect_equal(colSums(cat.SBS1536), colSums(x1))
  expect_equal("SBS96Catalog", class(x1)[1])

  cat.SBS1536.counts.signature <-
    TransformCatalog(cat.SBS1536, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")
  expect_error(Collapse1536CatalogTo96(cat.SBS1536.counts.signature))
})

test_that("Collapse192CatalogTo96 function is working properly", {
  cat.SBS192 <-
    ReadCatalog("testdata/regress.cat.sbs.192.csv", ref.genome = "GRCh37",
                region = "genome", catalog.type = "counts")
  x1 <- Collapse192CatalogTo96(cat.SBS192)
  expect_equal(colSums(cat.SBS192), colSums(x1))
  expect_equal("SBS96Catalog", class(x1)[1])

  cat.SBS192.counts.signature <-
    TransformCatalog(cat.SBS192, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")
  expect_error(Collapse192CatalogTo96(cat.SBS192.counts.signature))
})

test_that("Collapse144CatalogTo78 function is working properly", {
  cat.DBS144 <-
    ReadCatalog("testdata/regress.cat.dbs.144.csv", ref.genome = "GRCh37",
                region = "genome", catalog.type = "counts")
  x1 <- Collapse144CatalogTo78(cat.DBS144)
  expect_equal(colSums(cat.DBS144), colSums(x1))
  expect_equal("DBS78Catalog", class(x1)[1])

  cat.DBS144.counts.signature <-
    TransformCatalog(cat.DBS144, target.ref.genome = "GRCh37",
                     target.region = "genome",
                     target.catalog.type = "counts.signature")
  expect_error(Collapse144CatalogTo78(cat.DBS144.counts.signature))
})
