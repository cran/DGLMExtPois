#'Summarizing COM-Poisson Fits
#'
#'These functions are all methods for class \code{"glm_CMP"} or
#'\code{summary.glm_CMP} objects.
#'
#'@param object an object of class \code{"glm_CMP"}, usually, a result of a call
#'  to \code{glm.CMP}.
#'@param x an object of class \code{"summary.glm_CMP"}, usually, a result of a
#'  call to \code{summary.glm_CMP}.
#'@inheritParams stats::print.summary.glm
#' @examples
#' Bids$size.sq <- Bids$size^2
#' fit <- glm.CMP(formula.mu = numbids ~ leglrest + rearest + finrest +
#'                whtknght + bidprem + insthold + size + size.sq + regulatn,
#'                formula.nu = numbids ~ 1, data = Bids)
#' summary(fit)
#'@name summary.glm_CMP
NULL

#' @rdname summary.glm_CMP
#' @export
summary.glm_CMP <- function (object, ...) {
  # matrix of beta coefficients
  coef.p <- unname(object$betas)
  s.err  <- se_betas_cmp(object)
  if (is.null(s.err)) {
    coef.table <- matrix(coef.p, ncol = 1)
    rownames(coef.table) <- names(object$betas)
    colnames(coef.table) <- c("Estimate")
  } else {
    tvalue <- coef.p / s.err
    pvalue <- 2 * stats::pnorm(abs(tvalue), lower.tail = FALSE)
    coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    rownames(coef.table) <- names(object$betas)
    colnames(coef.table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  }
  # matrix of delta coefficients
  coef.p <- unname(object$deltas)
  s.err  <- se_deltas_cmp(object)
  if (is.null(s.err)) {
    coef.table.d <- matrix(coef.p, ncol = 1)
    rownames(coef.table.d) <- names(object$deltas)
    colnames(coef.table.d) <- c("Estimate")
  } else {
    tvalue <- coef.p / s.err
    pvalue <- 2 * stats::pnorm(abs(tvalue), lower.tail = FALSE)
    coef.table.d <- cbind(coef.p, s.err, tvalue, pvalue)
    rownames(coef.table.d) <- names(object$deltas)
    colnames(coef.table.d) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  }
  keep <- match(c("call", "aic"),
                names(object),
                0L
  )
  result <- object[keep]
  others <- list(
    coeff_betas = coef.table,
    coeff_deltas = coef.table.d
  )
  result <- c(result, others)
  class(result) <- "summary.glm_CMP"
  result
}

se_betas_cmp <- function(object) {
  if (! is.null(object$matrix.mu)) {
    matrix.mu <- object$matrix.mu
  } else if (! is.null(object$model.mu)) {
    matrix.mu <- stats::model.matrix(attr(object$model.mu, "terms"),
                                     object$model.mu)
  } else {
    formula.mu <- stats::as.formula(object$call[["formula.mu"]])
    a.mu <- stats::model.frame(formula.mu, data = object$data)
    matrix.mu <- stats::model.matrix(attr(a.mu, "terms"), a.mu)
  }

  In_beta  <- t(matrix.mu) %*% diag(as.vector((object$fitted.values ^ 2) /
                  variances_cmp(object$lambdas, object$nus, object$maxiter_series, object$tol))) %*%
                  matrix.mu

  tryCatch(sqrt(diag(solve(In_beta))),
           error = function(e) NULL
  )
}

se_deltas_cmp <- function(object) {
  deltas <- object$deltas
  lambda <- object$lambdas
  nu <- object$nus

  if (! is.null(object$y)) {
    y <- object$y
  } else if (! is.null(object$model.mu)) {
    y <- stats::model.response(object$model.mu)
  } else {
    formula.mu <- stats::as.formula(object$call[["formula.mu"]])
    a.mu <- stats::model.frame(formula.mu, data = object$data)
    y <- stats::model.response(a.mu)
  }

  if (! is.null(object$matrix.nu)) {
    matrix.nu <- object$matrix.nu
  } else if (! is.null(object$model.nu)) {
    matrix.nu <- stats::model.matrix(attr(object$model.nu, "terms"),
                                     object$model.nu)
  } else {
    formula.nu <- stats::as.formula(object$call[["formula.nu"]])
    a.nu <- stats::model.frame(formula.nu, data = object$data)
    matrix.nu <- stats::model.matrix(attr(a.nu, "terms"), a.nu)
  }

  W1 <- diag(as.vector(nu^2*variances_lfact(lambda, nu, object$maxiter_series, object$tol)))
  W2 <- diag(as.vector(nu^2*
                         covars_lfact_y(lambda, nu, object$maxiter_series, object$tol)^2/variances_cmp(lambda, nu, object$maxiter_series, object$tol)))

  In_delta  <- t(matrix.nu) %*% (W1-W2) %*% matrix.nu

  Z <- tryCatch(sqrt(diag(solve(In_delta))),
           error = function(e) NULL
  )

  if (! is.null(Z))
    return(Z)

  loglik_red <- function(deltas) {
    nus <- exp(matrix.nu %*% deltas)
    - sum(object$weights * (y * log(object$lambdas) - nus * lfactorial(y) -
                              log(Z(object$lambdas, nus, object$maxiter_series, object$tol))))
  }

  if (length(deltas) > 1) {
    fit_red <- stats::nlm(loglik_red, deltas, hessian = TRUE)
  } else {
    fit_red <- stats::optim(deltas, loglik_red, hessian = TRUE, method = "BFGS")
  }

  tryCatch(sqrt(diag(solve(fit_red$hessian))),
           error = function(e) NULL
  )
}

#' @rdname summary.glm_CMP
#' @export
print.summary.glm_CMP <- function (x, digits = max(3, getOption("digits") - 3),
                                  signif.stars = getOption("show.signif.stars"),
                                  ...) {
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n",
      sep = ""
  )
  cat("Mean model coefficients (with log link):\n")
  if (ncol(x$coeff_betas) > 1) {
    stats::printCoefmat(x$coeff_betas,
                        digits = digits,
                        signif.stars = signif.stars,
                        na.print = "NA",
                        ...
    )
  } else {
    print(x$coeff_betas)
  }
  cat("\nDispersion model coefficients (with logit link):\n")
  if (ncol(x$coeff_deltas) > 1) {
    stats::printCoefmat(x$coeff_deltas,
                        digits = digits,
                        signif.stars = signif.stars,
                        na.print = "NA",
                        ...
    )
  } else {
    print(x$coeff_deltas)
  }
  cat("\n")
  cat("AIC:", format(signif(x$aic, digits)), "\n")
  invisible(x)
}
