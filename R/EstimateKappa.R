# Estimate kappa for von Mises from robust measures of dispersion
#' @export
#' @importFrom stats complete.cases integrate quantile splinefun

EstimateKappa <- function(X, stat = c("CLMS", "CMAD", "CLTS"), results_grid = NULL) {
  stat <- match.arg(stat)

  if(!is.circular(X)){stop("argument X must be entered as a vector of class circular")}

  X  <- conversion.circular(X, units = "radians", modulo="2pi", zero=0, rotation="counter")


  # Compute robust statistic
  t_val <- switch(stat,
                  CMAD = CircularMAD(X, h = FALSE),
                  CLMS = CircularLMS(X),
                  CLTS = CircularLTS(X))

  # Generate grid if missing
  if (is.null(results_grid)) results_grid <- generate_vonmises_grid()

  # Interpolation
  valid_idx <- complete.cases(results_grid[[stat]], results_grid$CSD)
  interp_fun <- splinefun(
    x = results_grid[[stat]][valid_idx][order(results_grid[[stat]][valid_idx])],
    y = results_grid$CSD[valid_idx][order(results_grid[[stat]][valid_idx])],
    method = "monoH.FC"
  )

  est_CSD <- interp_fun(t_val)
  est_kappa <- circular::A1inv(1 - est_CSD / 2)

  list(CSD = est_CSD, kappa = est_kappa)
}
