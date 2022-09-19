test_that("Maximum number of evaluations", {
  opt <- list(maxeval = 10)
  Bids$size.sq <- Bids$size ^ 2
  fit <- glm.hP(formula.mu = numbids ~ leglrest + rearest + finrest +
                  whtknght + bidprem + insthold + size + size.sq + regulatn,
                formula.gamma = numbids ~ 1, data = Bids, opts = opt)
  expect_true(fit$nloptr$iterations <= 10)
})
