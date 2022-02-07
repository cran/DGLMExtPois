stopping <- function (x, tol) {
  all(abs(x) <= tol, na.rm = TRUE)
}
# stopping <- function (x, tol) all(abs(x) <= tol) & all(!is.nan(x))

# In the following functions a tolerance value less than 0 can be used to
# avoid the stopping criterion based on tolerance

# 1f1(1,gamma;lambda) to normalize hP
f11 <- function(lambda, gamma, maxiter_series = 10000, tol = 1.0e-10) {
  fac  <- 1
  temp <- 1
  L    <- gamma
  for (n in seq_len(maxiter_series)) {
    fac    <- fac * lambda / L
    series <- temp + fac
    if (stopping(series - temp, tol)){
      return(Re(series))
    }
    temp   <- series
    L      <- L + 1
  }
  if (tol >= 0)
    warning("Tolerance is not met")
  return(Re(series))
}

# E[Y]
means_hp <- function(lambda, gamma, maxiter_series = 10000, tol = 1.0e-10, f11_comp = NULL) {
  if (is.null(f11_comp))
    f11_comp <- f11(lambda, gamma, maxiter_series, tol)
  L    <- gamma
  fac  <- 1
  temp <- 0
  for (n in seq_len(maxiter_series)) {
    fac    <- fac * lambda / L
    series <- temp + n * fac
    if (stopping(series - temp, tol)){
      return(Re(series) / f11_comp)
    }
    temp   <- series
    L      <- L + 1
  }
  if (tol >= 0)
    warning("Tolerance is not met")
  return(Re(series) / f11_comp)
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
means_psiy <- function(lambda, gamma, maxiter_series = 10000, tol = 1.0e-10,
                       f11_comp = NULL) {
  if (is.null(f11_comp))
    f11_comp <- f11(lambda, gamma, maxiter_series, tol)
  L    <- gamma
  fac  <- 1
  temp <- digamma(gamma)
  for (n in seq_len(maxiter_series)) {
    fac    <- fac * lambda / L
    series <- temp + digamma(gamma + n) * fac
    if (stopping(series - temp, tol)){
      return(Re(series) / f11_comp)
    }
    temp   <- series
    L      <- L + 1
  }
  if (tol >= 0)
    warning("Tolerance is not met")
  return(Re(series) / f11_comp)
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

# Var(psi(gamma + Y))
variances_psiy <- function(lambda, gamma, maxiter_series = 10000, tol = 1.0e-10)
{
  f11_comp <- f11(lambda, gamma, maxiter_series, tol)
  L    <- gamma
  fac  <- 1
  temp <- digamma(gamma) ^ 2
  for (n in seq_len(maxiter_series)) {
    fac    <- fac * lambda / L
    series <- temp + digamma(gamma + n) ^ 2 * fac
    if (stopping(series - temp, tol)){
      return(Re(series) / f11_comp - means_psiy(lambda, gamma, maxiter_series, tol, f11_comp) ^ 2)
    }
    temp   <- series
    L      <- L + 1
  }
  if (tol >= 0)
    warning("Tolerance is not met")
  Re(series) / f11_comp - means_psiy(lambda, gamma, maxiter_series, tol, f11_comp) ^ 2
}

# Normalizing CMP
Z <- function(lambda, nu, maxiter_series, tol) {
  fac  <- 1
  temp <- 1
  for (n in seq_len(maxiter_series)) {
    fac    <- fac * lambda / (n^nu)
    series <- temp + fac
    if (stopping(series - temp, tol)){# & n >= 100){
      return(Re(series))
    }
    temp   <- series
  }
  return(Re(series))
}

# E[Y]
means_cmp <- function(lambda, nu, maxiter_series = 10000, tol = 1.0e-10) {
  fac  <- 1
  temp <- 0
  for (n in seq_len(maxiter_series)) {
    fac    <- fac * lambda / (n ^ nu)
    series <- temp + n * fac
    temp   <- series
  }
  return(Re(series) / Z(lambda, nu, maxiter_series, tol))
}

# E[log(Y!)]
means_lfact <- function(lambda, nu, maxiter_series = 10000, tol = 1.0e-10) {
  fac  <- 1
  temp <- 0
  for (n in seq_len(maxiter_series)) {
    fac    <- fac * lambda / (n ^ nu)
    series <- temp + lfactorial(n) * fac
    temp   <- series
  }
  return(Re(series) / Z(lambda, nu, maxiter_series, tol))
}

# E[Y log(Y!)]
means_lfact_y <- function(lambda, nu, maxiter_series = 10000, tol = 1.0e-10) {
  fac  <- 1
  temp <- 0
  for (n in seq_len(maxiter_series)) {
    fac    <- fac * lambda / (n^nu)
    series <- temp + n * lfactorial(n) * fac
    temp   <- series
  }
  return(Re(series) / Z(lambda, nu, maxiter_series, tol))
}

# Var(Y)
variances_cmp <- function(lambda, nu, maxiter_series = 10000, tol = 1.0e-10) {
  fac  <- 1
  temp <- 0
  for (n in seq_len(maxiter_series)) {
    fac    <- fac * lambda / (n^nu)
    series <- temp + n^2 * fac
    temp   <- series
  }
  return(Re(series) / Z(lambda, nu, maxiter_series, tol) - means_cmp(lambda, nu, maxiter_series, tol)^2)
}

# Var(log(Y!))
variances_lfact <- function(lambda, nu, maxiter_series = 10000, tol = 1.0e-10)
{
  fac  <- 1
  temp <- 0
  for (n in seq_len(maxiter_series)) {
    fac    <- fac * lambda / (n^nu)
    series <- temp + lfactorial(n) ^ 2 * fac
    temp   <- series
  }
  return(Re(series) / Z(lambda, nu, maxiter_series, tol) - means_lfact(lambda, nu, maxiter_series, tol) ^ 2)
}

# Cov(Y, log(Y!))
covars_lfact_y <- function(lambda, nu, maxiter_series = 10000, tol = 1.0e-10) {
  means_lfact_y(lambda, nu, maxiter_series, tol) - means_cmp(lambda, nu, maxiter_series, tol) * means_lfact(lambda, nu, maxiter_series, tol)
}
