test_that("fastmerge() concatenates dataframes", {
  x <- data.frame("X" = 1:2, "Y" = c("a", "b"))
  y <- data.frame("X" = 3:4, "Y" = c("a", "b"))
  expect_equal(fastmerge(x, y), data.frame("X" = 1:4, "Y" = c("a", "b", "a", "b")))
})
