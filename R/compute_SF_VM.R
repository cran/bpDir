# Compute CLTS population value for von Mises
#' @keywords internal
#' @importFrom circular dvonmises
#' @noRd
compute_SF_VM <- function(kappa) {
  lower <- qvonmises(0.25, mu = circular(0), kappa = kappa)
  upper <- qvonmises(0.75, mu = circular(0), kappa = kappa)

  f <- function(u) (1 - cos(u)) * dvonmises(circular(u), mu = circular(0), kappa = kappa)
  I <- integrate(f, lower, upper)$value
  2 * sqrt(I)
}
