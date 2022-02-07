#'Summarizing hyper-Poisson Fits
#'
#'These functions are all methods for class \code{"glm_hP"} or
#'\code{summary.glm_hP} objects.
#'
#'@param object an object of class \code{"glm_hP"}, usually, a result of a call
#'  to \code{glm.hP}.
#'@param x an object of class \code{"summary.glm_hP"}, usually, a result of a
#'  call to \code{summary.glm_hP}.
#'@inheritParams stats::print.summary.glm
#' @examples
#' ## Fit a hyper-Poisson model
#'
#' Bids$size.sq <- Bids$size ^ 2
#' fit <- glm.hP(formula.mu = numbids ~ leglrest + rearest + finrest +
#'               whtknght + bidprem + insthold + size + size.sq + regulatn,
#'               formula.gamma = numbids ~ 1, data = Bids)
#'
#' ## Obtain a summary of the fitted model
#'
#' summary(fit)
#'@name summary.glm_hP
NULL

#' @rdname summary.glm_hP
#' @export
summary.glm_hP <- function (object, ...) {
  # matrix of beta coefficients
  coef.p <- unname(object$betas)
  s.err  <- se_betas_hp(object)
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
  s.err  <- se_deltas_hp(object)
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
  class(result) <- "summary.glm_hP"
  result
}

se_betas_hp <- function(object) {
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

  n <- nrow(matrix.mu)
  betas   <- object$betas
  lambdas <- object$lambdas
  gammas  <- object$gammas
  offset  <- ifelse(is.null(object$offset), rep.int(0, n), object$offset)
  mu      <- exp(offset + matrix.mu %*% betas)

  In_beta  <- t(matrix.mu) %*% diag(as.vector((mu^2) /
              variances_hp(lambdas, gammas, object$maxiter_series, object$tol))) %*% matrix.mu
  tryCatch(sqrt(diag(solve(In_beta))),
           error = function(e) NULL
  )
}

se_deltas_hp <- function(object) {
  deltas <- object$deltas
  lambdas <- object$lambdas
  gammas <- object$gammas

  if (! is.null(object$y)) {
    y <- object$y
  } else if (! is.null(object$model.mu)) {
    y <- stats::model.response(object$model.mu)
  } else {
    formula.mu <- stats::as.formula(object$call[["formula.mu"]])
    a.mu <- stats::model.frame(formula.mu, data = object$data)
    y <- stats::model.response(a.mu)
  }

  if (! is.null(object$matrix.gamma)) {
    matrix.gamma <- object$matrix.gamma
  } else if (! is.null(object$model.gamma)) {
    matrix.gamma <- stats::model.matrix(attr(object$model.gamma, "terms"),
                                        object$model.gamma)
  } else {
    formula.gamma <- stats::as.formula(object$call[["formula.gamma"]])
    a.gamma <- stats::model.frame(formula.gamma, data = object$data)
    matrix.gamma <- stats::model.matrix(attr(a.gamma, "terms"), a.gamma)
  }

  W1 <- diag(as.vector(gammas^2*variances_psiy(lambdas, gammas, object$maxiter_series, object$tol)))
  W2 <- diag(as.vector(gammas^2*
                         covars_psiy(lambdas, gammas, object$maxiter_series, object$tol)^2/variances_hp(lambdas, gammas, object$maxiter_series, object$tol)))

  In_delta  <- t(matrix.gamma) %*% (W1-W2) %*% matrix.gamma

  Z <- tryCatch(sqrt(diag(solve(In_delta))),
           error = function(e) NULL
  )

  if (! is.null(Z))
    return(Z)

  loglik_red <- function(deltas) {
  gammas  <- exp(matrix.gamma %*% deltas)
    - sum(object$weights * (y * log(lambdas) -
                              lgamma(gammas + y) + lgamma(gammas) -
                              log(f11(lambdas, gammas, object$maxiter_series, object$tol))))
  }

  if (length(deltas) > 1){
    fit_red <- stats::optim(deltas, loglik_red, hessian = TRUE)

  } else {
    fit_red <- stats::optim(deltas, loglik_red, hessian = TRUE, method = "BFGS")
    }
  tryCatch(sqrt(diag(solve(fit_red$hessian))),
           error = function(e) NULL
  )
}

#' @rdname summary.glm_hP
#' @export
print.summary.glm_hP <- function (x, digits = max(3, getOption("digits") - 3),
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


#' @export
vcov.glm_hP <- function(object, ...) {
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

  n <- nrow(matrix.mu)
  betas   <- object$betas
  lambdas <- object$lambdas
  gammas  <- object$gammas
  offset  <- ifelse(is.null(object$offset), rep.int(0, n), object$offset)
  mu      <- exp(offset + matrix.mu %*% betas)

  In_beta  <- t(matrix.mu) %*% diag(as.vector((mu^2) /
              variances_hp(lambdas, gammas, object$maxiter_series, object$tol))) %*% matrix.mu
  solve(In_beta)
}
