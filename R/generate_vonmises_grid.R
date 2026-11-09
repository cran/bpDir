# Generate grid for interpolation for von Mises
#' @export
#' @importFrom circular qvonmises
#' @noRd
generate_vonmises_grid <- function(rho_seq = seq(0.01, 0.99, by = 0.002)) {
  kappa_vals <- circular::A1inv(rho_seq)
  results <- data.frame(
    kappa = kappa_vals,
    rho   = rho_seq,
    CSD     = 2 * (1 - rho_seq),
    CMAD  = NA,
    CLMS  = NA,
    CLTS  = NA
  )

  for (i in seq_along(rho_seq)) {
    kappa <- kappa_vals[i]
    results$CMAD[i] <- qvonmises(0.75, mu = circular(0), kappa = kappa)
    results$CLMS[i] <- qvonmises(0.75, mu = circular(0), kappa = kappa)
    results$CLTS[i] <- compute_SF_VM(kappa)
  }
  results
}
