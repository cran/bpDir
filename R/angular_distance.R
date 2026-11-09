#' @keywords internal
#' @noRd
#'
# Shortest angular distance between two directions on the circle.
angular_distance <- function(a1, a2) {
  d <- abs(a1 - a2)
  pmin(d, 2*pi - d)
}
