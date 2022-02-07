test_that("hP density function", {
  expect_equal(dhP(3, 15, 2), 0.0017017682)
  expect_equal(dhP(rep(3, 2), 15, 2), rep(0.0017017682, 2))
})
