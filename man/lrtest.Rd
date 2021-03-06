% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lrtest.R
\name{lrtest}
\alias{lrtest}
\title{Likelihood Ratio Test for Nested glm_CMP and glm_hP Fits}
\usage{
lrtest(object1, object2)
}
\arguments{
\item{object1, object2}{fitted objects of classes inheriting from
\code{"glm_CMP"} or \code{"glm_hP"}}
}
\value{
A list with class \code{"lrtest"} containing the following components:

  \item{\code{statistics}}{the value of the statistic.} \item{\code{df}}{the
  degrees of freedom.} \item{\code{p-value}}{the p-value for the test.}
}
\description{
Performs the likelihood ratio chi-squared test to compare nested models.
}
\details{
The test statistics is calculated as \eqn{2(llik- llik_0)}{2*(llik- llik_0)}.
The test statistics has a chi-squared distribution with r degrees of freedom,
where r is the difference in the number of parameters between the full and
null models.
}
\examples{
Bids$size.sq <- Bids$size ^ 2

## Fit null model
fit0 <- glm.hP(formula.mu = numbids ~ leglrest + rearest + finrest +
              whtknght + bidprem + insthold + size + size.sq + regulatn,
              formula.gamma = numbids ~ 1, data = Bids)

## Fit full model
fit <- glm.hP(formula.mu = numbids ~ leglrest + rearest + finrest +
              whtknght + bidprem + insthold + size + size.sq + regulatn,
              formula.gamma = numbids ~ leglrest, data = Bids)

## Likelihood ratio test for the nested models
lrtest(fit,fit0)
}
