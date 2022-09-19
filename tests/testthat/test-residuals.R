test_that("Absence of envelope", {
  Bids$size.sq <- Bids$size ^ 2
  hP.fit <- glm.hP(formula.mu = numbids ~ leglrest + rearest + finrest +
                 whtknght + bidprem + insthold + size + size.sq + regulatn,
                 formula.gamma = numbids ~ 1, data = Bids)
  expect_equal(residuals(hP.fit, "pearson", FALSE), res_hp(hP.fit, "pearson"))
  Bids$size.sq <- Bids$size ^ 2
  CMP.fit <- glm.CMP(formula.mu = numbids ~ leglrest + rearest + finrest +
                 whtknght + bidprem + insthold + size + size.sq + regulatn,
                 formula.nu = numbids ~ 1, data = Bids)
  expect_equal(residuals(CMP.fit, "pearson", FALSE), res_cmp(CMP.fit, "pearson"))
})
