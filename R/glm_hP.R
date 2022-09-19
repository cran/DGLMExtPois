#' Fit a hyper-Poisson Double Generalized Linear Model
#'
#' The \code{glm.hP} function is used to fit a hyper-Poisson double generalized
#' linear model with a log-link for the mean (\code{mu}) and the dispersion
#' parameter (\code{gamma}).
#'
#' Fit a hyper-Poisson double generalized linear model using as optimizer the
#' NLOPT_LD_SLSQP algorithm of function \code{\link[nloptr]{nloptr}}.
#'
#' @param formula.mu an object of class "formula" (or one that can be coerced to
#'   that class): a symbolic description of the model to be fitted.
#' @param formula.gamma regression formula linked to \code{log(gamma)}
#' @param init.beta initial values for regression coefficients of \code{beta}.
#' @param init.delta initial values for regression coefficients of \code{delta}.
#' @param data an optional data frame, list or environment (or object that can
#'   be coerced by \code{\link[base]{as.data.frame}} to a data frame) containing
#'   the variables in the model. If not found in data, the variables are taken
#'   from \code{environment(formula)}, typically the environment from which
#'   \code{glm.hP} is called.
#' @param weights an optional vector of \sQuote{prior weights} to be used in the
#'   fitting process. Should be \code{NULL} or a numeric vector.
#' @param subset an optional vector specifying a subset of observations to be
#'   used in the fitting process.
#' @param na.action a function which indicates what should happen when the data
#'   contain \code{NAs}. The default is set by the \code{na.action} setting of
#'   \code{\link[base]{options}}, and is \code{\link[stats]{na.fail}} if that is
#'   unset. The \sQuote{factory-fresh} default is
#'   \code{\link[stats:na.fail]{na.omit}}. Another possible value is
#'   \code{NULL}, no action. Value \code{\link[stats:na.fail]{na.exclude}} can
#'   be useful.
#' @param maxiter_series Maximum number of iterations to perform in the
#'   calculation of the normalizing constant.
#' @param tol tolerance with default zero meaning to iterate until additional
#'   terms to not change the partial sum in the calculation of the normalizing
#'   constant.
#' @param offset this can be used to specify an a priori known component to be
#'   included in the linear predictor during fitting. This should be \code{NULL}
#'   or a numeric vector of length equal to the number of cases. One or more
#'   \code{\link[stats]{offset}} terms can be included in the formula instead or
#'   as well, and if more than one is specified their sum is used. See
#'   \code{\link[stats:model.extract]{model.offset}}.
#' @param opts a list with options to the optimizer,
#'   \code{\link[nloptr]{nloptr}}, that fits the model. See, the \code{opts}
#'   parameter of \code{\link[nloptr]{nloptr}} for further details.
#' @param model.mu a logical value indicating whether the \emph{mu model frame}
#'   should be included as a component of the returned value.
#' @param model.gamma a logical value indicating whether the \emph{gamma model
#'   frame} should be included as a component of the returned value.
#' @param x logical value indicating whether the mu model matrix used in the
#'   fitting process should be returned as a component of the returned value.
#' @param y logical value indicating whether the response vector used in the
#'   fitting process should be returned as a component of the returned value.
#' @param z logical value indicating whether the gamma model matrix used in the
#'   fitting process should be returned as a component of the returned value.
#'
#' @return \code{glm.hP} returns an object of class \code{"glm_hP"}. The
#'   function \code{\link[base]{summary}} can be used to obtain or print a
#'   summary of the results.
#'
#'   The generic accessor functions \code{\link[stats]{coef}},
#'   \code{\link[stats]{fitted.values}} and \code{\link[stats]{residuals}} can
#'   be used to extract various useful features of the value returned by
#'   \code{glm.hP}.
#'
#'   \code{weights} extracts a vector of weights, one for each case in the fit
#'   (after subsetting and \code{na.action}).
#'
#'   An object of class \code{"glm_hP"} is a list containing at least the
#'   following components:
#'
#'   \item{\code{coefficients}}{a named vector of coefficients.}
#'   \item{\code{residuals}}{the residuals, that is response minus fitted
#'   values.} \item{\code{fitted.values}}{the fitted mean values.}
#'   \item{\code{linear.predictors}}{the linear fit on link scale.}
#'   \item{\code{call}}{the matched call.} \item{\code{offset}}{the offset
#'   vector used.} \item{\code{weights}}{the weights initially supplied, a
#'   vector of \code{1s} if none were.} \item{\code{df.residual}}{the residual
#'   degrees of freedom.} \item{\code{df.null}}{the residual degrees of freedom
#'   for the null model.} \item{\code{y}}{if requested (the default) the y
#'   vector used.} \item{\code{matrix.mu}}{if requested, the mu model matrix.}
#'   \item{\code{matrix.gamma}}{if requested, the gamma model matrix.}
#'   \item{\code{model.mu}}{if requested (the default) the mu model frame.}
#'   \item{\code{model.gamma}}{if requested (the default) the gamma model
#'   frame.} \item{\code{nloptr}}{an object of class \code{"nloptr"} with the
#'   result returned by the optimizer \code{\link[nloptr]{nloptr}}}
#'
#' @export
#'
#' @references
#'
#' Antonio J. Saez-Castillo and Antonio Conde-Sanchez (2013). "A hyper-Poisson
#' regression model for overdispersed and underdispersed count data",
#' Computational Statistics & Data Analysis, 61, pp. 148--157.
#'
#' S. G. Johnson (2018). \href{https://CRAN.R-project.org/package=nloptr}{The
#' nlopt nonlinear-optimization package}
#'
#' @examples
#' ## Fit model
#' Bids$size.sq <- Bids$size ^ 2
#' fit <- glm.hP(formula.mu = numbids ~ leglrest + rearest + finrest +
#'               whtknght + bidprem + insthold + size + size.sq + regulatn,
#'               formula.gamma = numbids ~ 1, data = Bids)
#'
#' ## Summary of the model
#' summary(fit)
#'
#' ## To see the termination condition of the optimization process
#' fit$nloptr$message
#'
#' ## To see the number of iterations of the optimization process
#' fit$nloptr$iterations
glm.hP <- function(formula.mu, formula.gamma, init.beta = NULL,
                   init.delta = NULL, data, weights, subset, na.action,
                   maxiter_series = 1000, tol = 0, offset,
                   opts = NULL, model.mu = TRUE,
                   model.gamma = TRUE, x = FALSE, y = TRUE, z = FALSE) {

  stopifnot(is.logical(model.mu))
  stopifnot(is.logical(model.gamma))
  stopifnot(is.logical(x))
  stopifnot(is.logical(y))
  stopifnot(is.logical(z))
  ret_y <- y

  formula.mu = stats::as.formula(formula.mu)
  if (base::missing(data)) data <- environment(formula.mu)

  # Design matrix a.mu
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
    offset <- rep.int(0, length(y))
  }

  matrizmu <- stats::model.matrix(attr(a.mu, "terms"), a.mu)
  ncovars.loc <- ncol(matrizmu)

  # Design matrix a.gamma

  mf <- match.call(expand.dots = FALSE)
  m  <- match(c("formula.gamma", "data", "subset", "weights", "na.action",
                "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  names(mf)[match("formula.gamma", names(mf))] <- "formula"
  a.gamma <- eval.parent(mf)

  matrizgamma <- stats::model.matrix(attr(a.gamma, "terms"), a.gamma)
  ncovars.gamma <- ncol(matrizgamma)

  # Optimization
  if (is.null(init.delta)) init.delta <- rep.int(0, ncovars.gamma)

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

  weights <- as.vector(stats::model.weights(a.mu))
  if (is.null(weights)) {
    weights <- rep.int(1, length(y))
  } else {
    if (! is.numeric(weights))
      stop("'weights' must be a numeric vector")
    if (any(weights < 0))
      stop("negative weights not allowed")
  }

  lambda0 <- log(rep(sum(y * weights) / sum(weights), length(y)))# + stats::rnorm(length(y), 0, 0.1)
  param0  <- c(init.beta, lambda0, init.delta)

  # Input variables
  n  <- length(y)
  q1 <- ncol(matrizmu) # Number of coefficients of mean equation
  q2 <- ncol(matrizgamma) # Number of coefficients of dispersion equation

  global_lambda <- NULL
  global_gamma <- NULL
  f11_cache <- NULL
  # 1f1(1,gamma;lambda) to normalize hP
  f11 <- function(lambda, gamma, maxiter_series = 10000, tol = 1.0e-10) {
    if (identical(lambda, global_lambda) && identical(gamma, global_gamma))
      return(f11_cache)
    global_lambda <<- lambda
    global_gamma <<- gamma
    fac  <- 1
    temp <- 1
    L    <- gamma
    for (n in seq_len(maxiter_series)) {
      fac    <- fac * lambda / L
      series <- temp + fac
      if (stopping(series - temp, tol)){
        f11_cache <<- Re(series)
        return(f11_cache)
      }
      temp   <- series
      L      <- L + 1
    }
    if (tol >= 0)
      warning("Tolerance is not met")
    f11_cache <<- Re(series)
    return(f11_cache)
  }

  # E[Y]
  lambda_means_hp <- NULL
  gamma_means_hp <- NULL
  means_hp_cache <- NULL
  means_hp <- function(lambda, gamma, maxiter_series = 10000, tol = 1.0e-10, f11_comp = NULL) {
    if (identical(lambda, lambda_means_hp) && identical(gamma, gamma_means_hp))
      return(means_hp_cache)
    lambda_means_hp <<- lambda
    gamma_means_hp <<- gamma
    if (is.null(f11_comp))
      f11_comp <- f11(lambda, gamma, maxiter_series, tol)
    L    <- gamma
    fac  <- 1
    temp <- 0
    for (n in seq_len(maxiter_series)) {
      fac    <- fac * lambda / L
      series <- temp + n * fac
      if (stopping(series - temp, tol)){
        means_hp_cache <<- Re(series) / f11_comp
        return(means_hp_cache)
      }
      temp   <- series
      L      <- L + 1
    }
    if (tol >= 0)
      warning("Tolerance is not met")
    means_hp_cache <<- Re(series) / f11_comp
    return(means_hp_cache)
  }

  # Var(Y)
  variances_hp <- function(lambda, gamma, maxiter_series = 10000, tol = 1.0e-10,
                           f11_comp = NULL, means_hp_comp = NULL)
  {
    if (is.null(f11_comp))
      f11_comp <- f11(lambda, gamma, maxiter_series, tol)
    if (is.null(means_hp_comp))
      means_hp_comp <- means_hp(lambda, gamma, maxiter_series, tol, f11_comp)
    L    <- gamma
    fac  <- 1
    temp <- 0
    for (n in seq_len(maxiter_series)) {
      fac    <- fac * lambda / L
      series <- temp + n ^ 2 * fac
      if (stopping(series - temp, tol)){
        return(Re(series) / f11_comp - means_hp_comp ^ 2)
      }
      temp   <- series
      L      <- L + 1
    }
    if (tol >= 0)
      warning("Tolerance is not met")
    return(Re(series) / f11_comp - means_hp_comp ^ 2)
  }

  # E[psi(gamma + Y)]
  lambda_means_psiy <- NULL
  gamma_means_psiy <- NULL
  means_psiy_cache <- NULL
  means_psiy <- function(lambda, gamma, maxiter_series = 10000, tol = 1.0e-10,
                         f11_comp = NULL) {
    if (identical(lambda, lambda_means_psiy) && identical(gamma, gamma_means_psiy))
      return(means_psiy_cache)
    lambda_means_psiy <<- lambda
    gamma_means_psiy <<- gamma
    if (is.null(f11_comp))
      f11_comp <- f11(lambda, gamma, maxiter_series, tol)
    L    <- gamma
    fac  <- 1
    temp <- digamma(gamma)
    for (n in seq_len(maxiter_series)) {
      fac    <- fac * lambda / L
      series <- temp + digamma(gamma + n) * fac
      if (stopping(series - temp, tol)){
        means_psiy_cache <<- Re(series) / f11_comp
        return(means_psiy_cache)
      }
      temp   <- series
      L      <- L + 1
    }
    if (tol >= 0)
      warning("Tolerance is not met")
    means_psiy_cache <<- Re(series) / f11_comp
    return(means_psiy_cache)
  }

  # Cov(Y, psi(gamma + Y))
  covars_psiy <- function(lambda, gamma, maxiter_series = 10000, tol = 1.0e-10,
                          f11_comp = NULL, means_hp_comp = NULL) {
    if (is.null(f11_comp))
      f11_comp <- f11(lambda, gamma, maxiter_series, tol)
    if (is.null(means_hp_comp))
      means_hp_comp <- means_hp(lambda, gamma, maxiter_series, tol, f11_comp)
    L    <- gamma
    fac  <- 1
    temp <- 0
    for (n in seq_len(maxiter_series)) {
      fac    <- fac * lambda / L
      series <- temp + n * digamma(gamma + n) * fac
      if (stopping(series - temp, tol)){
        return(Re(series) / f11_comp -
                 means_hp_comp *
                 means_psiy(lambda, gamma, maxiter_series, tol, f11_comp))
      }
      temp   <- series
      L      <- L + 1
    }
    if (tol >= 0)
      warning("Tolerance is not met")
    return(Re(series) / f11_comp -
             means_hp_comp *
             means_psiy(lambda, gamma, maxiter_series, tol, f11_comp))
  }

  loglik <- function(param) {
    lambda     <- exp(param[(q1 + 1):(q1 + n)])
    beta_gamma <- param[(q1 + n + 1):(q1 + n + q2)]
    gamma      <- as.vector(exp(matrizgamma %*% beta_gamma))
    value <- - sum(weights * (y * log(lambda) -
                              lgamma(gamma + y) + lgamma(gamma) -
                              log(f11(lambda, gamma, maxiter_series = maxiter_series, tol = tol))))
    return(value)
  }

  # Gradient of the objective function
  loglik_grad <- function(param) {
    lambda     <- exp(param[(q1+1):(q1+n)])
    beta_gamma <- param[(q1+n+1):(q1+n+q2)]
    gamma      <- as.vector(exp(matrizgamma %*% beta_gamma))
    f11_comp <- f11(lambda, gamma, maxiter_series = maxiter_series, tol = tol)
    value <- -c(rep(0, q1),
                weights * (y / lambda - means_hp(lambda, gamma, maxiter_series, tol, f11_comp) / lambda) * lambda,
                (weights * (-digamma(gamma + y) +
                              means_psiy(lambda, gamma, maxiter_series, tol, f11_comp))) %*%
                  (gamma * matrizgamma))
    #t(t(matrizgamma) %*% diag(gamma))
    value
  }

  # Constraints
  constraints <- function(param) {
    beta     <- param[1:q1]
    lambda   <- exp(param[(q1+1):(q1+n)])
    beta_gam <- param[(q1+n+1):(q1+n+q2)]
    gam      <- as.vector(exp(matrizgamma %*% beta_gam))
    exp(offset + matrizmu %*% beta) - means_hp(lambda, gam, maxiter_series, tol)
  }

  # Gradient constraints
  constraints_grad <- function(param) {
    beta_mu    <- param[1:q1]
    mu         <- exp(offset + matrizmu %*% beta_mu)
    lambda     <- exp(param[(q1+1):(q1+n)])
    beta_gamma <- param[(q1+n+1):(q1+n+q2)]
    gamma      <- as.vector(exp(matrizgamma %*% beta_gamma))
    # gradbeta   <- t(matrizmu) %*% diag(as.vector(mu))
    gradbeta   <- t(as.vector(mu) * matrizmu)
    f11_comp <- f11(lambda, gamma, maxiter_series = maxiter_series, tol = tol)
    means_hp_comp <- means_hp(lambda, gamma, maxiter_series, tol, f11_comp)
    gradlambda <- - diag(variances_hp(lambda, gamma, maxiter_series, tol, f11_comp, means_hp_comp) / lambda) * lambda
    v <- covars_psiy(lambda, gamma, maxiter_series, tol, f11_comp, means_hp_comp)
    # gradgamma <- diag(v * gamma) %*% matrizgamma
    gradgamma <- v * gamma * matrizgamma
    cbind(t(gradbeta), t(gradlambda), gradgamma)
  }

  # Optimization process
  my_local_opts <- list(algorithm = 'NLOPT_LD_SLSQP',
                        xtol_rel = 0.01
  )

  my_opts <- list(algorithm = 'NLOPT_LD_SLSQP',
                  maxeval = 1000, # 100000,
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

  fit <- nloptr::nloptr(param0,
                        eval_f = loglik,
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
    logver = fit$objective,
    bic = 2 * (fit$objective) + (q1 + q2) * log(sum(weights)),
    call = match.call(),
    formula.mu = formula.mu,
    formula.gamma = formula.gamma,
    df.residual = sum(weights) - (q1 + q2),
    df.null =  sum(weights) - 2,
    lambdas = exp(fit$pars[(q1 + 1):(q1 + n)]),
    gammas = exp(matrizgamma %*% fit$pars[(q1 + n + 1):(q1 + n + q2)]),
    coefficients2 = fit$pars,
    data = data,
    weights = stats::setNames(weights, seq(weights)),
    code = fit$status
  )
  results$fitted.values <- as.vector(exp(offset + matrizmu %*% fit$pars[1:q1]))
  names(results$fitted.values) <- seq(results$fitted.values)
  results$linear.predictors <- as.vector(offset + matrizmu %*% fit$pars[1:q1])
  names(results$linear.predictors) <- seq(results$linear.predictors)
  results$coefficients <-  fit$pars[1:q1]
  names(results$coefficients) <- colnames(matrizmu)
  results$betas <- results$coefficients
  results$deltas <- fit$pars[(q1 + n + 1):(q1 + n + q2)]
  names(results$deltas) <- colnames(matrizgamma)
  results$residuals <- stats::setNames(y - results$fitted.values, seq(y))
  results$maxiter_series <- maxiter_series
  results$tol <- tol

  if (ret_y) results$y <- y
  if (x) results$matrix.mu    <- matrizmu
  if (z) results$matrix.gamma <- matrizgamma
  if (model.mu) results$model.mu <- a.mu
  if (model.gamma) results$model.gamma <- a.gamma
  class(results) <- "glm_hP"
  results
}
