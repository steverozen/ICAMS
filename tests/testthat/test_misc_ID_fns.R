
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
