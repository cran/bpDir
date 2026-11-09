#' @keywords internal
#' @noRd
#'
# Build one or more polygon regions where `keep` is TRUE.
# This handles wrap-around by merging TRUE runs at the beginning and end.
make_region_polygons <- function(theta_plot, r_outer, r_inner, keep) {
  n <- length(keep)

  if (n == 0 || !any(keep)) {
    return(data.frame(x = numeric(0), y = numeric(0), group = integer(0)))
  }

  change_idx <- which(diff(keep) != 0)
  starts <- c(1, change_idx + 1)
  ends <- c(change_idx, n)

  runs <- Map(function(s, e) {
    list(start = s, end = e, value = keep[s])
  }, starts, ends)

  out <- list()
  k <- 1

  if (
    length(runs) > 1 &&
    isTRUE(runs[[1]]$value) &&
    isTRUE(runs[[length(runs)]]$value)
  ) {
    merged_idx <- c(
      runs[[length(runs)]]$start:runs[[length(runs)]]$end,
      runs[[1]]$start:runs[[1]]$end
    )

    if (length(merged_idx) >= 2) {
      out[[k]] <- data.frame(
        x = c(
          r_outer[merged_idx] * cos(theta_plot[merged_idx]),
          rev(r_inner[merged_idx] * cos(theta_plot[merged_idx]))
        ),
        y = c(
          r_outer[merged_idx] * sin(theta_plot[merged_idx]),
          rev(r_inner[merged_idx] * sin(theta_plot[merged_idx]))
        ),
        group = k
      )
      k <- k + 1
    }

    if (length(runs) > 2) {
      for (i in 2:(length(runs) - 1)) {
        if (isTRUE(runs[[i]]$value)) {
          idx <- runs[[i]]$start:runs[[i]]$end

          if (length(idx) >= 2) {
            out[[k]] <- data.frame(
              x = c(
                r_outer[idx] * cos(theta_plot[idx]),
                rev(r_inner[idx] * cos(theta_plot[idx]))
              ),
              y = c(
                r_outer[idx] * sin(theta_plot[idx]),
                rev(r_inner[idx] * sin(theta_plot[idx]))
              ),
              group = k
            )
            k <- k + 1
          }
        }
      }
    }

  } else {
    for (run in runs) {
      if (isTRUE(run$value)) {
        idx <- run$start:run$end

        if (length(idx) >= 2) {
          out[[k]] <- data.frame(
            x = c(
              r_outer[idx] * cos(theta_plot[idx]),
              rev(r_inner[idx] * cos(theta_plot[idx]))
            ),
            y = c(
              r_outer[idx] * sin(theta_plot[idx]),
              rev(r_inner[idx] * sin(theta_plot[idx]))
            ),
            group = k
          )
          k <- k + 1
        }
      }
    }
  }

  if (length(out) == 0) {
    data.frame(x = numeric(0), y = numeric(0), group = integer(0))
  } else {
    bind_rows(out)
  }
}

