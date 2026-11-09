#' @keywords internal
#' @noRd
#'
# Convert input values to circular radians with modulo 2pi.
as_circular_rad <- function(x) {
  circular(as.numeric(x), units = "radians", modulo = "2pi")
}
