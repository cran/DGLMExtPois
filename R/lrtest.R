#' Likelihood Ratio Test for Nested glm_CMP and glm_hP Fits
#'
#' Performs the likelihood ratio chi-squared test to compare nested models.
#'
#' The test statistics is calculated as \eqn{2(llik- llik_0)}{2*(llik- llik_0)}.
#' The test statistics has a chi-squared distribution with r degrees of freedom,
#' where r is the difference in the number of parameters between the full and
#' null models.
#'
#' @param object1,object2 fitted objects of classes inheriting from
#'   \code{"glm_CMP"} or \code{"glm_hP"}
#'
#' @return A list with class \code{"lrtest"} containing the following components:
#'
#'   \item{\code{statistics}}{the value of the statistic.} \item{\code{df}}{the
#'   degrees of freedom.} \item{\code{p-value}}{the p-value for the test.}
#'
#' @export
#'
#' @examples
#' Bids$size.sq <- Bids$size ^ 2
#'
#' ## Fit null model
#' fit0 <- glm.hP(formula.mu = numbids ~ leglrest + rearest + finrest +
#'               whtknght + bidprem + insthold + size + size.sq + regulatn,
#'               formula.gamma = numbids ~ 1, data = Bids)
#'
#' ## Fit full model
#' fit <- glm.hP(formula.mu = numbids ~ leglrest + rearest + finrest +
#'               whtknght + bidprem + insthold + size + size.sq + regulatn,
#'               formula.gamma = numbids ~ leglrest, data = Bids)
#'
#' ## Likelihood ratio test for the nested models
#' lrtest(fit,fit0)

lrtest <- function(object1, object2) {
  if(!(inherits(object1, "glm_CMP") && inherits(object2, "glm_CMP") ||
       inherits(object1, "glm_hP") && inherits(object2, "glm_hP")))
    stop("The objects must be of class glm_CMP or glm_hP")
  df <- abs(object1$df.residual - object2$df.residual)
  lrt <- abs(2*(object1$logver - object2$logver))
  result <- list(statistic = lrt,
                 df = df,
                 p.value = stats::pchisq(lrt, df, lower.tail = FALSE)
  )
  class(result) <- "lrtest"
  result
}

#' @export
print.lrtest <- function(x, ...) {
  cat("Statistics:", x$statistic, "\n")
  cat("Degrees of freedom:", x$df, "\n")
  cat("p-value:", x$p.value, "\n")
}
