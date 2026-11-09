# Circular Least Median Spread (CLMS)
# Computes shortest arc containing 50% of the data
#' @export
CircularLMS <- function(X) {
  n <- length(X)
  h <- ceiling((n + 1) / 2)

  # Order the data points
  Ordtheta <- as.vector(X[order(X)])

  # Extend vector to compute half-sample distances
  Newtheta <- c(Ordtheta, (Ordtheta[1:(h-1)] + 2*pi))

  # Compute half-sample minimal distance
  CLMS <- 0.5 * min(diff(Newtheta, lag = h-1))
  CLMS
}
