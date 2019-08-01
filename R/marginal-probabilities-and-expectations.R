#'Expected Probabilities and Frequencies for the hyper-Poisson and COM-Poisson
#'Model
#'
#'The \code{hP_expected} and \code{CMP_expected} functions calculate the
#'probability distribution of the count response variable Y for each observation
#'and obtain the corresponding expected frequencies. It is an informal way to
#'evaluate the fit of the hP or CMP model by comparing the predicted
#'distribution of counts with the observed distribution.
#'
#'The average expected probabilities are computed as \deqn{\bar(Pr)(y=k) =
#'\frac{1}{n} \sum_{i=1}^n \widehat{Pr}(y_i = k | x_i)}
#'
#'The expected frequencies are obtained by multiplying by n.
#'
#'Two measures are offered for summarizing the comparison between expected and
#'observed frequencies: the sum of the absolute value of differences and the sum
#'of the square of differences (similar to the Pearson statistic of goodness of
#'fit).
#'
#'@param object a fitted object of class inheriting from \code{"glm_hP"} or
#'  \code{"glm_CMP"}.
#'
#'@return A list containing the following components:
#'
#'  \item{\code{frequencies}}{the expected counts for the hP or CMP fit.}
#'  \item{\code{observed_freq}}{the observed distribution.}
#'  \item{\code{probabilities}}{the expected distribution for the hP or CMP
#'  fit.} \item{\code{dif}}{sum of the absolute value of differences between
#'  \code{frequencies} and \code{observed_freq}.} \item{\code{chi2}}{sum of the
#'  square of differences between \code{frequencies} and \code{observed_freq}.}
#'
#'@references Hilbe, J. M. (2011).Negative Binomial Regression. (2nd ed.).
#'  Cambridge University Press.#'
#'
#'  Long, J. S. & Freese, J. (2014). Regression Models for Categorical Dependent
#'  Variables using STATA. (3rd ed.). Stata Press.
#'
#'@name expected
NULL

#' @rdname expected
#' @examples
#' Bids$size.sq <- Bids$size ^ 2
#' hP.fit <- glm.hP(formula.mu = numbids ~ leglrest + rearest + finrest +
#'                  whtknght + bidprem + insthold + size + size.sq + regulatn,
#'                  formula.gamma = numbids ~ 1, data = Bids)
#' hP_expected(hP.fit)
#' @export
hP_expected <- function(object) {
  stopifnot(inherits(object, "glm_hP"))

  gamma  <- object$gamma
  lambda <- object$lambda

  if (! is.null(object$y)) {
    y <- object$y
  } else if (! is.null(object$model.mu)) {
    y <- stats::model.response(object$model.mu)
  } else {
    formula.mu <- stats::as.formula(object$call[["formula.mu"]])
    a.mu <- stats::model.frame(formula.mu, data = object$data)
    y <- stats::model.response(a.mu)
  }

  row <- max(y) + 1
  col <- length(lambda)
  i <- 0:(row - 1)
  f0 <- matrix(0, col)
  prob <- matrix(0, nrow = row + 1, ncol = col)
  dist <- matrix(0, nrow = row + 1, ncol = col)
  for (j in seq_len(col)) {
    prob[1:row, j] <- dhP(i, gamma[j], lambda[j])
    prob[row + 1, j] <- phP(row - 1, gamma[j], lambda[j], lower.tail = FALSE)
    dist[, j] <- object$weights[j] * prob[, j]
  }

  estimated_freq <- rowSums(dist)
  estimated_prob <- estimated_freq / sum(estimated_freq)

  observed_freq <- table(factor(y, levels = 0:row))
  observed_prob <- observed_freq / sum(observed_freq)

  dif <- sum(abs(observed_prob - estimated_prob))

  chi2 <- sum((observed_freq[1:(row-1)] - estimated_freq[1:(row-1)]) ^ 2 / estimated_freq[1:(row-1)]) +
    (sum(observed_freq[row]) - sum(estimated_freq[row:(row+1)])) ^ 2 / sum(estimated_freq[row:(row+1)])
  list(
    frequencies = estimated_freq,
    observed_freq = observed_freq,
    probabilities = estimated_prob,
    dif = dif,
    chi2 = chi2
  )
}

#' @rdname expected
#' @examples
#' Bids$size.sq <- Bids$size ^ 2
#' CMP.fit <- glm.CMP(formula.mu = numbids ~ leglrest + rearest + finrest +
#'                    whtknght + bidprem + insthold + size + size.sq + regulatn,
#'                    formula.nu = numbids ~ 1, data = Bids)
#' CMP_expected(CMP.fit)
#' @export
CMP_expected <- function(object) {
  stopifnot(inherits(object, "glm_CMP"))

  lambda <- object$lambdas
  nu     <- object$nus

  if (! is.null(object$y)) {
    y <- object$y
  } else if (! is.null(object$model.mu)) {
    y <- stats::model.response(object$model.mu)
  } else {
    formula.mu <- stats::as.formula(object$call[["formula.mu"]])
    a.mu <- stats::model.frame(formula.mu, data = object$data)
    y <- stats::model.response(a.mu)
  }

  row <- max(y) + 1
  col <- length(lambda)
  i <- 0:(row - 1)
  f0 <- matrix(0, col)
  prob <- matrix(0, nrow = row + 1, ncol = col)
  dist <- matrix(0, nrow = row + 1, ncol = col)
  for (j in seq_len(col)) {
    prob[1:row, j] <- compoisson::dcom(i, lambda[j], nu[j])
    prob[row + 1, j] <-  1 - sum(prob[, j])
    if (sum(prob[1:row, j]) > 1){
      prob[1:row, j] <- prob[1:row, j]/sum(prob[1:row, j])
      prob[row + 1, j] <- 0
    }
    dist[, j] <- object$weights[j] * prob[, j]
  }

  estimated_freq <- rowSums(dist)
  estimated_prob <- estimated_freq / sum(estimated_freq)

  observed_freq <- table(factor(y, levels = 0:row))
  observed_prob <- observed_freq / sum(observed_freq)

  dif <- sum(abs(observed_prob - estimated_prob))

  chi2 <- sum((observed_freq[1:(row-1)] - estimated_freq[1:(row-1)]) ^ 2 / estimated_freq[1:(row-1)]) +
    (sum(observed_freq[row]) - sum(estimated_freq[row:(row+1)])) ^ 2 / sum(estimated_freq[row:(row+1)])
  list(
    frequencies = estimated_freq,
    observed_freq = observed_freq,
    probabilities = estimated_prob,
    dif = dif,
    chi2 = chi2
  )
}
