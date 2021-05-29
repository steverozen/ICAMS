unlink("Rplots.pdf")
test_that("Check sigprofiler-formatted indel catalog rownames",
          { 
            spIDCat <- 
              read.table("testdata/SigPro-Cat/Strelka.ID.GRCh37.s1.ID83.tsv",
                         header = TRUE, row.names = 1)
            expect_equal(catalog.row.headers.sp$ID83,rownames(spIDCat))
            unlink("Rplots.pdf")
          })



test_that("TransRownames.ID.PCAWG.SigPro function", 
          {
            inputRownames <- ICAMS::catalog.row.order$ID
            outputRownames <- TransRownames.ID.PCAWG.SigPro(inputRownames)
            expect_true(setequal(outputRownames,catalog.row.headers.sp$ID83))
          })


test_that("TransRownames.ID.SigPro.PCAWG function", 
          {
            inputRownames <- catalog.row.headers.sp$ID83
            outputRownames <- TransRownames.ID.SigPro.PCAWG(inputRownames)
            expect_true(setequal(outputRownames,ICAMS::catalog.row.order$ID))
          })
