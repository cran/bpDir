# Estimate sigma from a robust measure of circular dispersion
#' @export
#' @importFrom stats complete.cases integrate quantile splinefun

EstimateSigma <- function(X, stat = c("CLMS", "CMAD", "CLTS"), results_grid = NULL) {
  stat <- match.arg(stat)

  if(!is.circular(X)){stop("argument X must be entered as a vector of class circular")}

  X  <- conversion.circular(X, units = "radians", modulo="2pi", zero=0, rotation="counter")


  # Compute robust measures of dispersion from observed data
  t_val <- switch(stat,
                  CMAD = CircularMAD(X, h = FALSE),
                  CLMS = CircularLMS(X),
                  CLTS = CircularLTS(X))

  # Generate grid if missing
  if (is.null(results_grid)) results_grid <- generate_results_grid_sigma()

  # Interpolation
  valid_idx <- which(!is.na(results_grid[[stat]]))
  interp_fun <- splinefun(
    x = results_grid[[stat]][valid_idx][order(results_grid[[stat]][valid_idx])],
    y = results_grid$CSD[valid_idx][order(results_grid[[stat]][valid_idx])],
    method = "monoH.FC"
  )

  est_CSD <- interp_fun(t_val)
  est_sigma <- sqrt(-2 * log(1 - est_CSD^2 / 2))

  list(CSD = est_CSD, sigma = est_sigma)
}
