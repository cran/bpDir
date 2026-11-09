# Generate results grid for interpolation
#' @export
#' @importFrom circular qwrappednormal

generate_results_grid_sigma <- function(sigma_seq = seq(0.1, 1.5, by = 0.0005)) {
  results <- data.frame(
    sigma = sigma_seq,
    rho   = exp(-sigma_seq^2 / 2),
    CSD     = sqrt(2 * (1 - exp(-sigma_seq^2 / 2))),
    CMAD  = NA,
    CLMS  = NA,
    CLTS  = NA
  )

  for (i in seq_along(sigma_seq)) {
    rho <- results$rho[i]
    results$CMAD[i] <- qwrappednormal(0.75, mu = circular(0), rho = rho)
    results$CLMS[i] <- qwrappednormal(0.75, mu = circular(0), rho = rho)
    results$CLTS[i] <- compute_SF_WN(rho)
  }
  results
}
