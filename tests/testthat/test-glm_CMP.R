test_that("Default maximum iterations", {
  Bids$size.sq <- Bids$size^2
  fit <- glm.CMP(formula.mu = numbids ~ leglrest + rearest + finrest +
                  whtknght + bidprem + insthold + size + size.sq + regulatn,
                  formula.nu = numbids ~ 1, data = Bids)
  expect_equal(fit$maxiter_series, 1000)
})
