#' @keywords internal
#' @noRd
#'
# Build an annular arc polygon between start_angle and end_angle.
# Used for black box arcs, median patch, and CI interval.
arc_ribbon_df <- function(start_angle, end_angle, r_inner, r_outer, n = 200) {
  if (end_angle < start_angle) end_angle <- end_angle + 2*pi
  angles <- seq(start_angle, end_angle, length.out = n) %% (2*pi)

  data.frame(
    x = c(r_outer * cos(angles), rev(r_inner * cos(angles))),
    y = c(r_outer * sin(angles), rev(r_inner * sin(angles)))
  )
}
