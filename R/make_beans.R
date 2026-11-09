#' @keywords internal
#' @noRd
#' @importFrom dplyr bind_rows
#'
# Build bean segments at observation angles.
# Bean length is radial; linewidth is scaled by local violin thickness.
make_beans <- function(angles, theta_plot, r_outer, r_inner, side, bean_min_size, bean_max_size, bean_length) {
  if (length(angles) == 0) {
    return(data.frame(
      x = numeric(0), y = numeric(0),
      xend = numeric(0), yend = numeric(0),
      size = numeric(0)
    ))
  }

  local_thickness <- abs(r_outer - r_inner)

  # Keep only positions where the violin region has positive thickness
  valid <- which(local_thickness > 0)

  # Safety guard: if no valid positions exist, return empty output
  if (length(valid) == 0) {
    return(data.frame(
      x = numeric(0), y = numeric(0),
      xend = numeric(0), yend = numeric(0),
      size = numeric(0)
    ))
  }

  tmin <- min(local_thickness[valid], na.rm = TRUE)
  tmax <- max(local_thickness[valid], na.rm = TRUE)

  out <- lapply(seq_along(angles), function(i) {
    th <- angles[i]

    # Find nearest plotting angle only among valid positions
    j <- valid[
      which.min(angular_distance(theta_plot[valid], th))
    ]

    d_scaled <- (local_thickness[j] - tmin) / (tmax - tmin + 1e-9)
    size <- bean_min_size + d_scaled * (bean_max_size - bean_min_size)

    if (side == "inner") {
      r_start <- r_inner[j]
      r_end <- r_inner[j] - bean_length
    } else {
      r_start <- r_outer[j]
      r_end <- r_outer[j] + bean_length
    }

    data.frame(
      x = r_start * cos(th),
      y = r_start * sin(th),
      xend = r_end * cos(th),
      yend = r_end * sin(th),
      size = size
    )
  })

  dplyr::bind_rows(out)
}

