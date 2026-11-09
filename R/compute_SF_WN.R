# Compute CLTS poulation value for wrapped normal distribution
#' @keywords internal
#' @importFrom circular dwrappednormal qwrappednormal
#' @noRd
compute_SF_WN <- function(rho) {
  lower <- qwrappednormal(0.25, mu = circular(0), rho = rho)
  upper <- qwrappednormal(0.75, mu = circular(0), rho = rho)

  f <- function(u) (1 - cos(u)) * dwrappednormal(circular(u), mu = circular(0), rho = rho)
  I <- integrate(f, lower, upper)$value
  2 * sqrt(I)
}
