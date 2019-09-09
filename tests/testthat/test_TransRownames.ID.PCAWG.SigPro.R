

test_that("Check sigprofiler-formatted indel catalog rownames",
          {
            spIDCat <- read.csv("testdata/sigProfiler_ID_signatures.csv",
                            row.names = 1, 
                            stringsAsFactors = F)
            expect_equal(ICAMS::catalog.row.order.sp$ID83,rownames(spIDCat))
          })



#test_that("TransRownames.ID.PCAWG.SigPro function", 
#          {
#            inputRownames <- ICAMS::catalog.row.order$ID
#            outputRownames <- TransRownames.ID.PCAWG.SigPro(inputRownames)
#            expect_equal(outputRownames,ICAMS::catalog.row.order.sp$ID83)
#          })


#test_that("TransRownames.ID.SigPro.PCAWG function", 
#          {
#            inputRownames <- ICAMS::catalog.row.order.sp$ID83
#            outputRownames <- TransRownames.ID.PCAWG.SigPro(inputRownames)
#            expect_equal(outputRownames,ICAMS::catalog.row.order$ID)
#          })
