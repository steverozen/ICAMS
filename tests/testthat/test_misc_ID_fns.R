
context("Categorizing the mutation type of deletions")

test_that("FindMaxRepeatDel", {
  expect_equal(FindMaxRepeatDel("GATAGATAGATA", rep.unit.seq = "GATA", pos = 1), 2)
  expect_equal(FindMaxRepeatDel("GATAGATAGATA", rep.unit.seq = "GATA", pos = 5), 2)
  expect_equal(FindMaxRepeatDel("GATAGATAGATA", rep.unit.seq = "GATA", pos = 9), 2)
  expect_equal(FindMaxRepeatDel("CACAGATATATA", rep.unit.seq = "GATA", pos = 5), 0)
  expect_equal(FindMaxRepeatDel("CCCCCC", "C", pos = 3), 5)
  expect_error(FindMaxRepeatDel("CCCCCC", "T", pos = 3),
               "rep.unit.seq is not TRUE")
})

test_that("FindMaxRepeatIns", {
  tmp <- "xyaczt"
  expect_equal(FindMaxRepeatIns(tmp, "g", pos = 1), 0)
  expect_equal(FindMaxRepeatIns(tmp, "g", pos = 0), 0)
  expect_equal(FindMaxRepeatIns(tmp, "x", pos = 0), 1)
  expect_equal(FindMaxRepeatIns(tmp, "t", pos = 6), 1)
  expect_equal(FindMaxRepeatIns(tmp, "ac", pos = 2), 1)
  expect_equal(FindMaxRepeatIns(tmp, "ac", pos = 4), 1)
  for (i in c(1, 3, 5, 7, 9)) {
    expect_equal(FindMaxRepeatIns("gacacacacg", "ac", pos = i), 4)
  }
  for (i in c(0, 2, 4, 6, 8)) {
    expect_equal(FindMaxRepeatIns("acacacac", "ac", pos = i), 4)
  }
})
