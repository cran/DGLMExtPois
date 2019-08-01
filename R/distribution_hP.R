#' The hyper-Poisson Distribution
#'
#' Density, distribution function and random generation for the hyper-Poisson
#' distribution with parameters \code{gamma} and \code{lambda}.
#'
#' @param x vector of (non-negative integer) quantiles.
#' @param q vector of quantiles.
#' @param n	number of random values to return.
#' @param gamma dispersion parameter. Must be strictly positive.
#' @param lambda location parameter. Must be strictly positive.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return \code{dhP} gives the density, \code{phP} gives the distribution
#'   function and \code{rhP} generates random deviates.
#'
#'   Invalid \code{gamma} or \code{lambda} will result in return value
#'   \code{NaN}, with a warning.
#'
#'   The length of the result is determined by n for \code{rhP}, and is the
#'   maximum of the lengths of the numerical arguments for the other functions.
#'
#'   The numerical arguments other than \code{n} are recycled to the length of
#'   the result. Only the first element of the logical arguments is used.
#' @name hP
NULL

#' @rdname hP
#' @export
#'
#' @examples
#' dhP(3, 15, 2)
dhP <-function(x, gamma, lambda) {
  if (!(is.double(x) || is.integer(x)) ||
      !(is.double(gamma) || is.integer(gamma)) ||
      !(is.double(lambda) || is.integer(lambda)))
    stop("Non-numeric argument to mathematical function")
  maximum <- max(length(x), length(gamma), length(lambda))
  x <- rep(x, length.out = maximum)
  gamma  <- rep(gamma, length.out = maximum)
  lambda <- rep(lambda, length.out = maximum)
  fj <- numeric(length = maximum)
  warn <- FALSE
  for (ind in seq_along(x)) {
    if (gamma[ind] <= 0 || lambda[ind] <= 0) {
      fj[ind] <- NaN
      warn <- TRUE
    } else if (!(abs(x[ind]-round(x[ind]))<.Machine$double.eps)) {
      warning(paste("non-integer x =", x[ind]))
      fj[ind] <- 0
    }
    else if (x[ind] < 0) {
      fj[ind] <- 0
    } else if (x[ind] == 0) {
      fj[ind] <- 1 / f11(lambda[ind], gamma[ind])
    } else if (x[ind] == 1) {
      f0 <- 1 / f11(lambda[ind], gamma[ind])
      fj[ind] <- f0 * (f11(lambda[ind],gamma[ind],x[ind], tol = -1)-1)
    }
    else {
      f0 <- 1 / f11(lambda[ind], gamma[ind])
      fj[ind] <- f0 * (f11(lambda[ind],gamma[ind],x[ind], tol = -1) -
                         f11(lambda[ind],gamma[ind],x[ind]-1, tol = -1))
    }
  }
  if (warn)
    warning("NaN(s) produced: gamma and lambda must be strictly positive")
  return(fj)
}

#' @rdname hP
#' @export
#'
#' @examples
#' phP(3, 15, 2)
phP <- function(q, gamma, lambda, lower.tail = TRUE) {
  if (!(is.double(q) || is.integer(q)) ||
      !(is.double(gamma) || is.integer(gamma)) ||
      !(is.double(lambda) || is.integer(lambda)))
    stop("Non-numeric argument to mathematical function")
  maximum <- max(length(q), length(gamma), length(lambda))
  q <- rep(q, length.out = maximum)
  gamma  <- rep(gamma, length.out = maximum)
  lambda <- rep(lambda, length.out = maximum)
  prob <- numeric(length = maximum)
  for (ind in seq_along(q)) {
    if (q[ind] < 0) {
      prob[ind] <- 0
      next
    }
    qlower <- trunc(q[ind])
    x <- 0:qlower
    probs <- dhP(x, gamma[ind], lambda[ind])
    if (lower.tail) {
      prob[ind] <- sum(probs)
    } else {
      prob[ind] <- 1 - sum(probs)
    }
  }
  return(prob)
}

#' @rdname hP
#' @export
#'
#' @examples
#' rhP(10, 15, 2)
rhP <-function(n, gamma, lambda) {
  # check n parameter
  if(!is.numeric(n) || length(n) != 1 || n < 0)
    stop("invalid arguments")
  if (!(is.double(gamma) || is.integer(gamma)) ||
      !(is.double(lambda) || is.integer(lambda)))
    stop("Non-numeric argument to mathematical function")
  gamma  <- rep(gamma, length.out = n)
  lambda <- rep(lambda, length.out = n)
  result <- numeric(length = n)
  warn <- FALSE
  for (ind in seq_len(n)) {
    if (gamma[ind] <= 0 || lambda[ind] <= 0) {
      result[ind] <- NaN
      warn <- TRUE
    } else {
      result[ind] <- simulate_hp(gamma[ind], lambda[ind])
    }
  }
  if (warn)
    warning("NaN(s) produced: gamma and lambda must be strictly positive")
  result
}
