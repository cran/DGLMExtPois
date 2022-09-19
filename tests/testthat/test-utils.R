test_that("stopping function works", {
  expect_equal(stopping(c(0.2, 4), 1), FALSE)
  expect_equal(stopping(c(0.2, 0.5), 1), TRUE)
})
