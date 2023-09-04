#'AIC and BIC for COM-Poisson Fitted Models
#'
#'Computes the Akaike's information criterion or the Schwarz's Bayesian criterion
#'for COM-Poisson Fits
#'
#'@param object an object of class \code{"glm_CMP"}, typically the result of a
#'  call to \code{\link{glm.CMP}}.
#'@inheritParams stats::AIC
#'@name AIC_CMP
NULL

#' @rdname AIC_CMP
#' @export
#' @examples
#' ## Estimate a COM-Poisson model
#' Bids$size.sq <- Bids$size ^ 2
#' fit <- glm.CMP(formula.mu = numbids ~ leglrest + rearest + finrest +
#'               whtknght + bidprem + insthold + size + size.sq + regulatn,
#'               formula.nu = numbids ~ 1, data = Bids)
#'
#' ## Compute its AIC and BIC
#' AIC(fit)
#' @export
AIC.glm_CMP <- function(object, ..., k = 2) {
  others <- list(...)
  if (length(others) == 0)
    return(object$aic)
  lapply(others, function(x) stopifnot(inherits(x, "glm_CMP")))
  c(object$aic, sapply(others, function(x) x$aic))
}

#' @rdname AIC_CMP
#' @export
#' @examples
#' BIC(fit)
#' @export
BIC.glm_CMP <- function(object, ...) {
  others <- list(...)
  if (length(others) == 0)
    return(object$bic)
  lapply(others, function(x) stopifnot(inherits(x, "glm_CMP")))
  c(object$aic, sapply(others, function(x) x$bic))
}

#' @export
print.glm_CMP <- function (x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:  ",
      paste(deparse(x$call),
            sep = "\n",
            collapse = "\n"
      ),
      "\n\n",
      sep = ""
  )
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits),
                print.gap = 2, quote = FALSE)
  cat("\nDegrees of Freedom:",
      x$df.null,
      "Total (i.e. Null);",
      x$df.residual,
      " Residual\n"
  )
  cat("AIC:", format(signif(x$aic, digits))
  )
  invisible(x)
}

#' Predict Method for glm_CMP Fits
#'
#' Obtains predictions from a fitted \code{glm_CMP} object.
#'
#' @param object	a fitted object of class inheriting from \code{"glm_CMP"}.
#' @param newdata	optionally, a data frame in which to look for variables with
#'   which to predict. If omitted, the fitted linear predictors are used.
#' @param type the type of prediction required. The default is on the scale of
#'   the linear predictors; the alternative \code{"response"} is on the scale of
#'   the response variable.
#' @param ... 	further arguments passed to or from other methods.
#'
#' @return A vector with the prediction means.
#'
#' @examples
#' ## Fit a model
#' Bids$size.sq <- Bids$size ^ 2
#' fit <- glm.CMP(formula.mu = numbids ~ leglrest + rearest + finrest +
#'                whtknght + bidprem + insthold + size + size.sq + regulatn,
#'                formula.nu = numbids ~ 1, data = Bids)
#'
#' ## As the newdata parameter is not used the fitted values are obtained
#' predict(fit, type = "response")
#' @export
predict.glm_CMP <- function(object, newdata = NULL,
                           type = c("link", "response"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    r <- switch (type,
                 "link" = object$linear.predictors,
                 "response" = object$fitted.values,
                 stop("Invalid type value")
    )
    return(r)
  }
  formula.mu <- stats::as.formula(object$call[["formula.mu"]])
  a.mu <- stats::model.frame(formula.mu, data = newdata)
  offset <- stats::model.extract(a.mu, "offset")
  if (is.null(offset)) {
    offset <- rep.int(0, length(nrow(newdata)))
  }
  matrizmu <- stats::model.matrix(stats::terms(formula.mu),
                                  stats::model.frame(stats::terms(formula.mu),
                                                     data = newdata,
                                                     na.action = NULL)
  )
  if (type == "link") {
    r <- offset + matrizmu %*% object$todo$pars[1:ncol(matrizmu)]
  } else if (type == "response") {
    r <- exp(offset + matrizmu %*% object$todo$pars[1:ncol(matrizmu)])
  }
  r <- as.vector(r)
  names(r) <- seq(r)
  r
}

#' @export
coef.glm_CMP <- function(object, ...) {
  list(mean_model = object$betas,
       dispersion_model = object$deltas
  )
}

#' Confidence Intervals for glm_hP Fits
#'
#' Computes confidence intervals for one or more parameters in a \code{glm_CMP}
#' object.
#'
#' @param object a fitted object of class inheriting from \code{"glm_CMP"}.
#' @param parm a specification of which parameters are to be given confidence
#'   intervals, either a vector of numbers or a vector of names. If missing, all
#'   parameters are considered.
#' @param level the confidence level required.
#' @param ... additional argument(s) for methods.
#' @return A matrix (or vector) with columns giving lower and upper confidence
#'   limits for each \code{beta} parameter. These will be labelled as
#'   (1-level)/2 and 1 - (1-level)/2 in (by default 2.5\% and 97.5\%).
#' @examples
#' ## Estimate the model
#' Bids$size.sq <- Bids$size ^ 2
#' fit <- glm.CMP(formula.mu = numbids ~ leglrest + rearest + finrest +
#'                whtknght + bidprem + insthold + size + size.sq + regulatn,
#'                formula.nu = numbids ~ 1, data = Bids)
#'
#' ## Compute confidence intervals
#' confint(fit)
#' @export
confint.glm_CMP <- function (object, parm, level = 0.95, ...)
{
  cf <- object$betas
  pnames <- names(cf)
  if (missing(parm)) parm <- pnames
  else if(is.numeric(parm)) parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- format_perc(a, 3)
  fac <- stats::qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L),
              dimnames = list(parm, pct))
  ses <- se_betas_cmp(object)[parm]
  ci[] <- cf[parm] + ses %o% fac
  ci
}
