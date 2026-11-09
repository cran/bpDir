#' @keywords internal
#' @noRd
#' @importFrom ggplot2 geom_segment
#'
# Add one geometry object to the ggplot.
# This is called once for full/outer/inner and twice for two_halves.
add_geom <- function(
    p,
    geom,
    fill,
    alpha,
    fill_outlier,
    alpha_outlier,
    arc_color,
    ci_color,
    beans,
    bean_min_size,
    bean_max_size,
    bean_length,
    notch,
    show_notch_cutout,
    show_notch_guides,
    r0,
    outlier_point_color,
    outlier_point_size,
    outlier_point_alpha
) {

  # Draw the outlier density region for every beans option except
  # "outlier_points". In that mode, outliers are shown as points instead
  # of as an outlier density region.
  if (beans != "outlier_points" && nrow(geom$df_poly_out) > 0) {
    p <- p +
      geom_polygon(
        data = geom$df_poly_out,
        aes(x, y, group = group),
        fill = fill_outlier,
        color = "black",
        alpha = alpha_outlier,
        linewidth = 0.35
      )
  }

  if (nrow(geom$df_poly_in) > 0) {
    p <- p +
      geom_polygon(
        data = geom$df_poly_in,
        aes(x, y, group = group),
        fill = fill,
        color = "black",
        alpha = alpha,
        linewidth = 0.4
      )
  }

  bean_side <- ifelse(geom$side == "inner", "inner", "outer")

  if (beans %in% c("all", "inliers")) {
    df_beans <- make_beans(
      geom$inlier_angles_plot,
      geom$theta_plot,
      geom$r_outer_full,
      geom$r_inner_full,
      bean_side,
      bean_min_size,
      bean_max_size,
      bean_length
    )

    if (nrow(df_beans) > 0) {
      p <- p +
        geom_segment(
          data = df_beans,
          aes(x = x, y = y, xend = xend, yend = yend, linewidth = size),
          color = "black",
          alpha = 0.8,
          lineend = "butt"
        )
    }
  }

  if (beans %in% c("all", "outliers")) {
    df_out_beans <- make_beans(
      geom$outlier_angles_plot,
      geom$theta_plot,
      geom$r_outer_out,
      geom$r_inner_out,
      bean_side,
      bean_min_size,
      bean_max_size,
      bean_length
    )

    if (nrow(df_out_beans) > 0) {
      p <- p +
        geom_segment(
          data = df_out_beans,
          aes(x = x, y = y, xend = xend, yend = yend, linewidth = size),
          color = "black",
          alpha = 0.65,
          lineend = "butt"
        )
    }
  }

  # Draw the notch outline before the black quartile box.
  # This prevents the notch boundary from sitting on top of the box edge.
  # The black box is drawn afterward and masks the notch termination cleanly.
  if (notch && show_notch_cutout && nrow(geom$df_notch_cutout) > 0) {
    p <- p +
      geom_polygon(
        data = geom$df_notch_cutout,
        aes(x, y, group = group),
        fill = NA,
        color = "black",
        linewidth = 0.4
      )
  }

  # Draw the black quartile box arcs after the notch outline.
  p <- p +
    geom_polygon(data = geom$df_arc_Q1, aes(x, y), fill = arc_color) +
    geom_polygon(data = geom$df_arc_Q3, aes(x, y), fill = arc_color)

  # Draw the median patch after the quartile arcs.
  p <- p +
    geom_polygon(
      data = geom$df_median,
      aes(x, y),
      fill = fill,
      color = NA
    )

  if (notch && show_notch_guides && nrow(geom$df_notch_guides) > 0) {
    p <- p +
      geom_path(
        data = geom$df_notch_guides,
        aes(x, y, group = group),
        color = "black",
        linewidth = notch_guide_linewidth,
        linetype = notch_guide_linetype
      )
  }

  if (
    beans == "outlier_points" &&
    length(geom$outlier_angles_plot) > 0
  ) {
    df_out_points <- data.frame(
      x = r0 * cos(geom$outlier_angles_plot),
      y = r0 * sin(geom$outlier_angles_plot)
    )

    p <- p +
      geom_point(
        data = df_out_points,
        aes(x, y),
        color = outlier_point_color,
        size = outlier_point_size,
        alpha = outlier_point_alpha
      )
  }

  p
}

