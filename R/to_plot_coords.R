#' @keywords internal
#' @noRd
#'
# Convert mathematical angles to plotting angles using the chosen zero and rotation.
to_plot_coords <- function(theta, rotation, zero) {
  if (rotation == "counter") {
    (theta + zero) %% (2*pi)
  } else {
    (zero - theta) %% (2*pi)
  }
}
