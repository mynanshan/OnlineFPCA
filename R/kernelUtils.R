kernel.rbf <- function(d, h) {
  return(exp(-d^2 / h))
}

kernel.matern <- function(d, nu, phi, sigma2 = 1, .eps = 1e-10) {
  out <- d
  if (nu == 1) stop("NOT IMPLEMENTED")
  # nonzero part
  x <- out[abs(d) >= .eps]
  z <- sqrt(2 * nu) * x / phi
  out[abs(d) >= .eps] <-
    sigma2 * exp(log((z)^nu * besselK(z, nu)) - log(2^(nu - 1) * gamma(nu)))
  # limits at zero
  x <- out[abs(d) < .eps]
  z <- x / phi
  out[abs(d) < .eps] <-
    sigma2 * (1 + nu / (2 * (1 - nu)) * z^2 + nu^2 / (8 * (2 - 3 * nu + nu^2)) * z^4)
  return(out)
}
