#' @keywords internal
#' @noRd
#'
# Build the notch outline polygon.
# For outer: pmax(r_notched, notch_limit) clips the notch at the outer box boundary.
# For inner: pmin(r_notched, notch_limit) clips the notch at the inner box boundary.
# The V shape comes from r_notched, which was tapered around the median.
make_notch_cutout_df <- function(theta_plot, idx, r_full, r_notched, notch_limit, side) {
  if (length(idx) < 2) {
    return(data.frame(x = numeric(0), y = numeric(0), group = integer(0)))
  }

  if (side == "outer") {
    keep <- r_full[idx] > notch_limit
    if (sum(keep) < 2) {
      return(data.frame(x = numeric(0), y = numeric(0), group = integer(0)))
    }

    idx2 <- idx[keep]
    r_outer_use <- r_full[idx2]
    r_inner_use <- pmax(r_notched[idx2], notch_limit)

  } else if (side == "inner") {
    keep <- r_full[idx] < notch_limit
    if (sum(keep) < 2) {
      return(data.frame(x = numeric(0), y = numeric(0), group = integer(0)))
    }

    idx2 <- idx[keep]
    r_outer_use <- pmin(r_notched[idx2], notch_limit)
    r_inner_use <- r_full[idx2]

  } else {
    return(data.frame(x = numeric(0), y = numeric(0), group = integer(0)))
  }

  data.frame(
    x = c(
      r_outer_use * cos(theta_plot[idx2]),
      rev(r_inner_use * cos(theta_plot[idx2]))
    ),
    y = c(
      r_outer_use * sin(theta_plot[idx2]),
      rev(r_inner_use * sin(theta_plot[idx2]))
    ),
    group = 1
  )
}
