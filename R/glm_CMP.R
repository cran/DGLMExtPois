#' Fit a COM-Poisson Double Generalized Linear Model
#'
#' The \code{glm.CMP} function is used to fit a COM-Poisson double generalized
#' linear model with a log-link for the mean (\code{mu}) and the dispersion
#' parameter (\code{nu}).
#'
#' Fit a COM-Poisson double generalized linear model using as optimizer the
#' NLOPT_LD_SLSQP algorithm of function \code{\link[nloptr]{nloptr}}.
#'
#' @param formula.mu regression formula linked to \code{log(mu)}
#' @param formula.nu regression formula linked to \code{log(nu)}
#' @param init.beta initial values for regression coefficients of \code{beta}.
#' @param init.delta initial values for regression coefficients of \code{delta}.
#' @param data an optional data frame, list or environment (or object coercible
#'   by \code{\link[base]{as.data.frame}} to a data frame) containing the
#'   variables in the model. If not found in data, the variables are taken from
#'   \code{environment(formula)}, typically the environment from which
#'   \code{glm.CMP} is called.
#' @param model.nu a logical value indicating whether the \emph{nu model frame}
#'   should be included as a component of the returned value.
#' @param z logical value indicating whether the nu model matrix used in the
#'   fitting process should be returned as a component of the returned value.
#' @inheritParams glm.hP
#'
#' @return \code{glm.CMP} returns an object of class \code{"glm_CMP"}. The
#'   function \code{\link[base]{summary}} can be used to obtain or print a
#'   summary of the results. An object of class \code{"glm_CMP"} is a list
#'   containing at least the following components: \item{\code{coefficients}}{a
#'   named vector of coefficients.} \item{\code{residuals}}{the residuals, that
#'   is response minus fitted values.} \item{\code{fitted.values}}{the fitted
#'   mean values.} \item{\code{linear.predictors}}{the linear fit on link
#'   scale.} \item{\code{call}}{the matched call.} \item{\code{offset}}{the
#'   offset vector used.} \item{\code{weights}}{the weights initially supplied,
#'   a vector of \code{1s} if none were.} \item{\code{y}}{if requested (the
#'   default) the y vector used.} \item{\code{matrix.mu}}{if requested, the mu
#'   model matrix.} \item{\code{matrix.nu}}{if requested, the nu model matrix.}
#'   \item{\code{model.mu}}{if requested (the default) the mu model frame.}
#'   \item{\code{model.nu}}{if requested (the default) the nu model
#'   frame.}\item{\code{nloptr}}{an object of class \code{"nloptr"} with the
#'   result returned by the optimizer \code{\link[nloptr]{nloptr}}}
#' @export
#'
#' @references Alan Huang (2017). "Mean-parametrized Conway–Maxwell–Poisson
#'   regression models for dispersed counts", Statistical Modelling, 17(6), pp.
#'   359--380.
#'
#'   S. G. Johnson (2018). \href{https://CRAN.R-project.org/package=nloptr}{The
#'   nlopt nonlinear-optimization package}
#'
#' @examples
#' ## Fit model
#' Bids$size.sq <- Bids$size^2
#' fit <- glm.CMP(formula.mu = numbids ~ leglrest + rearest + finrest +
#'                whtknght + bidprem + insthold + size + size.sq + regulatn,
#'                formula.nu = numbids ~ 1, data = Bids)
#'
#' ## Summary of the model
#' summary(fit)
#'
#' ## To see termination condition of the optimization process
#' fit$nloptr$message
#'
#' ## To see number of iterations of the optimization process
#' fit$nloptr$iterations
glm.CMP <- function(formula.mu, formula.nu, init.beta = NULL,
                    init.delta = NULL, data, weights, subset, na.action,
                    maxiter_series = 1000, tol = 0, offset, opts = NULL,
                    model.mu = TRUE, model.nu = TRUE, x = FALSE, y = TRUE,
                    z = FALSE) {

  stopifnot(is.logical(model.mu))
  stopifnot(is.logical(model.nu))
  stopifnot(is.logical(x))
  stopifnot(is.logical(y))
  stopifnot(is.logical(z))
  ret_y <- y

  # Design matrix mu ------------------------------------------------------
  mf <- match.call(expand.dots = FALSE)
  m  <- match(c("formula.mu", "data", "subset", "weights", "na.action",
                "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  names(mf)[match("formula.mu", names(mf))] <- "formula"
  a.mu <- eval.parent(mf)
  offset   <- stats::model.extract(a.mu, "offset")
  y <- stats::model.extract(a.mu, "response")
  if (is.null(offset)) {
    offset    <- rep.int(0, length(y))
  }
  matrizmu <- stats::model.matrix(attr(a.mu, "terms"), a.mu)

  # Design matrix nu ------------------------------------------------------
  mf <- match.call(expand.dots = FALSE)
  m  <- match(c("formula.nu", "data", "subset", "weights", "na.action",
                "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  names(mf)[match("formula.nu", names(mf))] <- "formula"
  a.nu <- eval.parent(mf)
  matriznu <- stats::model.matrix(attr(a.nu, "terms"), a.nu)

  if (is.null(init.delta)) init.delta <- rep.int(0, ncol(matriznu))

  if (is.null(init.beta)) {
    parameters <- list(
      formula = formula.mu,
      data = data,
      family = "poisson"
    )
    if (! missing(subset)) parameters$subset <- subset
    if (! missing(na.action)) parameters$na.action <- na.action
    if (! missing(weights)) parameters$weights <- weights
    init.beta <- do.call(stats::glm, parameters)$coefficients
    }

  weights <- stats::model.extract(a.mu, "weights")
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  } else {
    if (! is.numeric(weights))
      stop("'weights' must be a numeric vector")
    if (any(weights < 0))
      stop("negative weights not allowed")
  }

  global_lambda_Z <- NULL
  global_nu_Z <- NULL
  Z_cache <- NULL
  Z <- function(lambda, nu, maxiter_series, tol) {
    if (identical(lambda, global_lambda_Z) && identical(nu, global_nu_Z))
      return(Z_cache)
    global_lambda_Z <<- lambda
    global_nu_Z <<- nu
    fac  <- 1
    temp <- 1
    for (n in seq_len(maxiter_series)) {
      fac    <- fac * lambda / (n^nu)
      series <- temp + fac
      if (stopping(series - temp, tol)){  # & n >= 100){
        Z_cache <<- Re(series)
        return(Z_cache)
      }
      temp <- series
    }
    Z_cache <<- Re(series)
    return(Z_cache)
  }

  global_lambda_means_cmp <- NULL
  global_nu_means_cmp <- NULL
  means_cmp_cache <- NULL
  means_cmp <- function(lambda, nu, maxiter_series = 10000, tol = 1.0e-10) {
    if (identical(lambda, global_lambda_means_cmp) && identical(nu, global_nu_means_cmp))
      return(means_cmp_cache)
    global_lambda_means_cmp <<- lambda
    global_nu_means_cmp <<- nu
    fac  <- 1
    temp <- 0
    for (n in seq_len(maxiter_series)) {
      fac    <- fac * lambda / (n ^ nu)
      series <- temp + n * fac
      temp   <- series
    }
    means_cmp_cache <<- Re(series) / Z(lambda, nu, maxiter_series, tol)
    means_cmp_cache
  }

  global_lambda_means_lfact <- NULL
  global_nu_means_lfact <- NULL
  means_lfact_cache <- NULL
  means_lfact <- function(lambda, nu, maxiter_series = 10000, tol = 1.0e-10) {
    if (identical(lambda, global_lambda_means_lfact) && identical(nu, global_nu_means_lfact))
      return(means_lfact_cache)
    global_lambda_means_lfact <<- lambda
    global_nu_means_lfact <<- nu
    fac  <- 1
    temp <- 0
    for (n in seq_len(maxiter_series)) {
      fac    <- fac * lambda / (n ^ nu)
      series <- temp + lfactorial(n) * fac
      temp   <- series
    }
    means_lfact_cache <<- Re(series) / Z(lambda, nu, maxiter_series, tol)
    means_lfact_cache
  }

  global_lambda_variances_cmp <- NULL
  global_nu_variances_cmp <- NULL
  variances_cmp_cache <- NULL
  variances_cmp <- function(lambda, nu, maxiter_series = 10000, tol = 1.0e-10) {
    if (identical(lambda, global_lambda_variances_cmp) &&
        identical(nu, global_nu_variances_cmp))
      return(variances_cmp_cache)
    global_lambda_variances_cmp <<- lambda
    global_nu_variances_cmp <<- nu
    fac  <- 1
    temp <- 0
    for (n in seq_len(maxiter_series)) {
      fac    <- fac * lambda / (n^nu)
      series <- temp + n^2 * fac
      temp   <- series
    }
    variances_cmp_cache <<- Re(series) / Z(lambda, nu, maxiter_series, tol) - means_cmp(lambda, nu, maxiter_series, tol)^2
    variances_cmp_cache
  }

  n  <- length(y)
  q1 <- ncol(matrizmu)
  q2 <- ncol(matriznu)

  loglik <- function(param) {
    lambda  <- exp(param[(q1+1): (q1+n)])
    beta_nu <- param[(q1+n+1):(q1+n+q2)]
    nu      <- as.vector(exp(matriznu %*% beta_nu))
    return(- sum(weights * (y * log(lambda) - nu * lfactorial(y) -
                              log(Z(lambda, nu, maxiter_series, tol)))))
  }

  loglik_grad <- function(param) {
    lambda  <- exp(param[(q1+1): (q1+n)])
    beta_nu <- param[(q1+n+1):(q1+n+q2)]
    nu      <- as.vector(exp(matriznu %*% beta_nu))
    return(-c(rep(0, q1),
              weights * (y / lambda - means_cmp(lambda, nu, maxiter_series, tol) / lambda) * lambda,
              (weights * (-lfactorial(y) + means_lfact(lambda, nu, maxiter_series, tol))) %*%
                t(t(matriznu) %*% diag(nu))))
  }

  constraints <- function(param) {
    beta    <- param[1:q1]
    lambda  <- exp(param[(q1+1):(q1+n)])
    beta_nu <- param[(q1+n+1):(q1+n+q2)]
    nu      <- as.vector(exp(matriznu %*% beta_nu))
    return(exp(offset + matrizmu %*% beta) - means_cmp(lambda, nu, maxiter_series, tol))
  }

  constraints_grad  <- function(param) {
    beta    <- param[1:q1]
    lambda  <- exp(param[(q1+1):(q1+n)])
    beta_nu <- param[(q1+n+1):(q1+n+q2)]
    nu      <- as.vector(exp(matriznu %*% beta_nu))
    gradbeta <- t(matrizmu) %*% diag(as.vector(exp(offset + matrizmu %*% beta)))
    gradlambda <- - diag(variances_cmp(lambda, nu, maxiter_series, tol) / lambda) * lambda
    gradnu     <- diag((means_lfact_y(lambda, nu, maxiter_series, tol) - means_cmp(lambda, nu, maxiter_series, tol) *
                          means_lfact(lambda, nu, maxiter_series, tol)) * nu) %*% matriznu
    return(cbind(t(gradbeta), t(gradlambda), gradnu))
  }

  # Optimization

  # starting values for optimization
  delta0  <- rep(log(1), q2)
  #lambda0 <- log(rep(mean(y), n))
  lambda0 <- log(rep(sum(y * weights) / sum(weights), length(y)))
  #lambda0 <- log(exp(offset + matrizmu %*% init.beta))
  param0  <- c(init.beta, lambda0, delta0)

  my_local_opts <- list(algorithm = 'NLOPT_LD_SLSQP',
                        xtol_rel = 0.01
  )
  my_opts <- list(algorithm = 'NLOPT_LD_SLSQP',
                  maxeval = 1000,
                  local_opts = my_local_opts,
                  print_level = 0
  )
  if (! is.null(opts)) {
    if ("local_opts" %in% names(opts)) {
      my_opts$local_opts[names(opts$local_opts)] <- opts$local_opts
      opts$local_opts <- NULL
    }
    my_opts[names(opts)] <- opts
  }
  fit <- nloptr::nloptr(param0, eval_f = loglik,
                        eval_grad_f = loglik_grad,
                        eval_g_eq = constraints,
                        eval_jac_g_eq = constraints_grad,
                        opts = my_opts
  )
  fit$pars <- fit$solution

  results <- list(
    nloptr = fit,
    offset = unname(stats::model.extract(a.mu, "offset")),
    aic = 2 * (fit$objective) + (q1 + q2) * 2,
    bic = 2 * (fit$objective) + (q1 + q2) * log(sum(weights)),
    logver = fit$objective,
    df.residual = sum(weights) - (q1 + q2),
    df.null =  sum(weights) - 2,
    call = match.call(),
    formula.mu = formula.mu,
    formula.nu = formula.nu,
    lambdas = exp(fit$pars[(q1 + 1):(q1 + n)]),
    nus = as.vector(exp(matriznu %*% fit$pars[(q1 + n + 1):(q1 + n + q2)])),
    coefficients2 = fit$pars,
    coefficients = stats::setNames(fit$pars[1:q1], colnames(matrizmu)),
    data = data,
    weights = stats::setNames(weights, seq(weights)),
    code = fit$status
  )

  results$fitted.values <- as.vector(exp(offset + matrizmu %*% fit$pars[1:q1]))
  names(results$fitted.values) <- seq(results$fitted.values)
  results$linear.predictors <- as.vector(offset + matrizmu %*% fit$pars[1:q1])
  names(results$linear.predictors) <- seq(results$linear.predictors)
  results$residuals <- stats::setNames(y - results$fitted.values, seq(y))
  results$betas <- fit$pars[1:q1]
  names(results$betas) <- colnames(matrizmu)
  results$deltas <- fit$pars[(q1+n+1):(q1+n+q2)]
  names(results$deltas) <- colnames(matriznu)
  results$maxiter_series <- maxiter_series
  results$tol <- tol

  if (ret_y) results$y <- y
  if (x) results$matrix.mu <- matrizmu
  if (z) results$matrix.nu <- matriznu
  if (model.mu) results$model.mu <- a.mu
  if (model.nu) results$model.nu <- a.nu
  class(results) <- "glm_CMP"
  results
}
