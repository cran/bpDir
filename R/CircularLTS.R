#' @export
CircularLTS <- function(X) {
  # Number of observations
  n <- length(X)
  h <- ceiling((n + 1) / 2)

  # Order the angles
  Ordtheta <- as.vector(X[order(X)])

  # Extend the ordered sample to handle circular wrap-around.
  Newtheta <- c(Ordtheta, Ordtheta[1:(h-1)] + 2*pi)

  # Store mean resultant lengths for each candidate half-sample
  Rhos <- numeric(n)

  # Slide a window of fixed size h over all consecutive subsets
  # Each subset corresponds to an arc containing h observations,
  for (i in 1:n) {
    subset <- Newtheta[i:(i + h - 1)]

    # Compute mean resultant length rho of this subset
    # Minimizing circular variance is equivalent to maximizing rho
    Rhos[i] <- rho.circular(circular(subset))
  }

  # The CLTS is the minimum CSD over all admissible subsets,
  # which is equivalent to using the maximum rho.
  sqrt(2 * (1 - max(Rhos)))
}
