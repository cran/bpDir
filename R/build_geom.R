#' @keywords internal
#' @noRd
#' @importFrom circular density.circular
#'
# Core geometry builder for one dataset and one side mode.
# It computes densities, quartiles, outliers, CI, notch, polygons, and angles.
build_geom <- function(
    x,
    n_boot,
    ci_level,
    kernel,
    density_scale,
    density_threshold,
    notch,
    notch_factor,
    r0,
    #rotation,
    #zero,
    arc_thickness,
    box_width_factor,
    ci_thickness,
    side = c("both", "outer", "inner"),
    bw_use_arg = bw,          # bandwidth used for this geometry
    cutoff_use_arg = cutoff     # cutoff used for this geometry
) {

  #side <- match.arg(side)

  if (length(x) < 3) {
    stop("`data` must contain at least 3 observations.")
  }

  # Standardize the input to circular radians.
  x <- as_circular_rad(x)

  # Open the circle at the anti-median to avoid discontinuity problems.
  med <- as.numeric(median.circular(x))
  anti_med <- (med + pi) %% (2*pi)

  x_open <- as_circular_rad((as.numeric(x) - anti_med) %% (2*pi))
  data_angles <- as.numeric(x_open)

  # Compute quartiles and median on the opened circular scale.
  Q1 <- as.numeric(quantile(x_open, probs = 0.25))
  Median <- as.numeric(median.circular(x_open))
  Q3 <- as.numeric(quantile(x_open, probs = 0.75))

  # Circular fences used to split inlier and outlier regions.
  lower_fence <- (Median - cutoff_use_arg + 2*pi) %% (2*pi)
  upper_fence <- (Median + cutoff_use_arg) %% (2*pi)

  inside <- if (lower_fence < upper_fence) {
    data_angles[data_angles >= lower_fence & data_angles <= upper_fence]
  } else {
    data_angles[data_angles >= lower_fence | data_angles <= upper_fence]
  }

  if (length(inside) == 0) inside <- data_angles

  outliers <- setdiff(data_angles, inside)
  data_inliers <- as_circular_rad(inside)
  n_in <- length(data_inliers)

  # Bootstrap the circular median using inlier observations only.
  boot_medians <- replicate(
    n_boot,
    as.numeric(median.circular(sample(data_inliers, n_in, replace = TRUE)))
  )

  # Convert bootstrap medians into angular deviations from the observed median.
  boot_dev <- angular_distance(boot_medians, Median)
  ci_radius <- as.numeric(
    quantile(boot_dev, 1 - (1 - ci_level) / 2, na.rm = TRUE)
  )

  # Angular limits of the notch/CI interval.
  notch_lower <- (Median - ci_radius) %% (2*pi)
  notch_upper <- (Median + ci_radius) %% (2*pi)

  bw_use <- if (is.null(bw_use_arg)) stats::bw.nrd(data_angles) else bw_use_arg
  if (!is.finite(bw_use) || length(bw_use) != 1 || bw_use <= 0) bw_use <- 0.3

  # Circular kernel density estimate on [0, 2*pi].
  dens <- density.circular(
    x_open,
    bw = bw_use,
    kernel = kernel,
    from = 0,
    to = 2*pi,
    n = 1000
  )

  theta <- as.numeric(dens$x)
  # Normalize density to the chosen visual radial thickness.
  d_vals <- dens$y / max(dens$y, na.rm = TRUE) * density_scale
  d_vals[d_vals < density_threshold] <- 0

  min_angle <- min(inside)
  max_angle <- max(inside)

  in_range <- if (min_angle < max_angle) {
    theta >= min_angle & theta <= max_angle
  } else {
    theta >= min_angle | theta <= max_angle
  }

  out_range <- !in_range

  in_notch <- if (notch_lower < notch_upper) {
    theta >= notch_lower & theta <= notch_upper
  } else {
    theta >= notch_lower | theta <= notch_upper
  }

  # Taper equals 1 everywhere unless a notch is requested.
  taper <- rep(1, length(theta))

  # Apply V-shaped taper inside the notch interval when notch = TRUE.
  if (notch && any(in_notch & in_range)) {
    dist_from_median <- angular_distance(theta[in_notch & in_range], Median)
    max_dist <- max(dist_from_median)

    if (max_dist > 0) {
      taper[in_notch & in_range] <-
        notch_factor * (dist_from_median / max_dist) + (1 - notch_factor)
    } else {
      taper[in_notch & in_range] <- 1 - notch_factor
    }
  }

  # Split density into inlier, outlier, and tapered inlier components.
  d_in <- ifelse(in_range, d_vals, 0)
  d_out <- ifelse(out_range, d_vals, 0)
  d_notched <- d_in * taper

  # Convert density values to inner/outer radii according to side mode.
  if (side == "both") {
    r_outer_full <- r0 + d_in
    r_inner_full <- r0 - d_in

    r_outer_notched <- r0 + d_notched
    r_inner_notched <- r0 - d_notched

    r_outer_out <- r0 + d_out
    r_inner_out <- r0 - d_out
  }

  if (side == "outer") {
    r_outer_full <- r0 + d_in
    r_inner_full <- rep(r0, length(theta))

    r_outer_notched <- r0 + d_notched
    r_inner_notched <- rep(r0, length(theta))

    r_outer_out <- r0 + d_out
    r_inner_out <- rep(r0, length(theta))
  }

  if (side == "inner") {
    r_outer_full <- rep(r0, length(theta))
    r_inner_full <- r0 - d_in

    r_outer_notched <- rep(r0, length(theta))
    r_inner_notched <- r0 - d_notched

    r_outer_out <- rep(r0, length(theta))
    r_inner_out <- r0 - d_out
  }

  # Close angles back to the original circular coordinate system.
  close_back <- function(a) (a + anti_med) %% (2*pi)

  #theta_plot <- to_plot_coords(close_back(theta))#, rotation,zero)
  #Q1_plot <- to_plot_coords(close_back(Q1))#, rotation,zero)
  #Median_plot <- to_plot_coords(close_back(Median))#, rotation,zero)
  #Q3_plot <- to_plot_coords(close_back(Q3))#, rotation,zero)
  #notch_lower_plot <- to_plot_coords(close_back(notch_lower))#, rotation,zero)
  #notch_upper_plot <- to_plot_coords(close_back(notch_upper))#, rotation,zero)
  #inlier_angles_plot <- to_plot_coords(close_back(inside))#, rotation,zero)
  #outlier_angles_plot <- to_plot_coords(close_back(outliers))#, rotation,zero)

  theta_plot <- close_back(theta)
  Q1_plot <- close_back(Q1)
  Median_plot <- close_back(Median)
  Q3_plot <- close_back(Q3)
  notch_lower_plot <- close_back(notch_lower)
  notch_upper_plot <- close_back(notch_upper)
  inlier_angles_plot <- close_back(inside)
  outlier_angles_plot <- close_back(outliers)

  # Polygon for the inlier density region.
  df_poly_in <- make_region_polygons(
    theta_plot,
    r_outer_full,
    r_inner_full,
    d_in > 0
  )

  # Polygon for the outlier density region.
  df_poly_out <- make_region_polygons(
    theta_plot,
    r_outer_out,
    r_inner_out,
    d_out > 0
  )

  # Final thickness of the black quartile box.
  box_arc_thickness <- arc_thickness * box_width_factor

  if (side == "outer") {
    arc_r_inner <- r0
    arc_r_outer <- r0 + box_arc_thickness
    median_r_inner <- r0
    median_r_outer <- r0 + box_arc_thickness

  } else if (side == "inner") {
    arc_r_inner <- r0 - box_arc_thickness
    arc_r_outer <- r0
    median_r_inner <- r0 - box_arc_thickness
    median_r_outer <- r0

  } else {
    arc_r_inner <- r0 - box_arc_thickness / 2
    arc_r_outer <- r0 + box_arc_thickness / 2
    median_r_inner <- r0 - box_arc_thickness / 2
    median_r_outer <- r0 + box_arc_thickness / 2
  }

  # Density grid positions belonging to the notch interval.
  notch_idx <- which(in_notch & in_range & d_in > 0)

  i1 <- which.min(angular_distance(theta_plot, notch_lower_plot))
  i2 <- which.min(angular_distance(theta_plot, notch_upper_plot))

  # Outer-side notch is clipped at arc_r_outer, the outer boundary of the black box.
  if (side == "outer") {
    df_notch_cutout <- make_notch_cutout_df(
      theta_plot = theta_plot,
      idx = notch_idx,
      r_full = r_outer_full,
      r_notched = r_outer_notched,
      notch_limit = arc_r_outer,
      side = "outer"
    )
    df_notch_cutout$group <- 1

    df_notch_guides <- data.frame(
      x = c(
        r_outer_notched[i1] * cos(notch_lower_plot),
        arc_r_outer * cos(notch_lower_plot),
        r_outer_notched[i2] * cos(notch_upper_plot),
        arc_r_outer * cos(notch_upper_plot)
      ),
      y = c(
        r_outer_notched[i1] * sin(notch_lower_plot),
        arc_r_outer * sin(notch_lower_plot),
        r_outer_notched[i2] * sin(notch_upper_plot),
        arc_r_outer * sin(notch_upper_plot)
      ),
      group = c(1, 1, 2, 2)
    )

    # Inner-side notch is clipped at arc_r_inner, the inner boundary of the black box.
  } else if (side == "inner") {
    df_notch_cutout <- make_notch_cutout_df(
      theta_plot = theta_plot,
      idx = notch_idx,
      r_full = r_inner_full,
      r_notched = r_inner_notched,
      notch_limit = arc_r_inner,
      side = "inner"
    )
    df_notch_cutout$group <- 1

    df_notch_guides <- data.frame(
      x = c(
        r_inner_notched[i1] * cos(notch_lower_plot),
        arc_r_inner * cos(notch_lower_plot),
        r_inner_notched[i2] * cos(notch_upper_plot),
        arc_r_inner * cos(notch_upper_plot)
      ),
      y = c(
        r_inner_notched[i1] * sin(notch_lower_plot),
        arc_r_inner * sin(notch_lower_plot),
        r_inner_notched[i2] * sin(notch_upper_plot),
        arc_r_inner * sin(notch_upper_plot)
      ),
      group = c(1, 1, 2, 2)
    )

    # Full violin uses both outer and inner notch outlines.
  } else {
    notch_outer <- make_notch_cutout_df(
      theta_plot = theta_plot,
      idx = notch_idx,
      r_full = r_outer_full,
      r_notched = r_outer_notched,
      notch_limit = arc_r_outer,
      side = "outer"
    )

    notch_inner <- make_notch_cutout_df(
      theta_plot = theta_plot,
      idx = notch_idx,
      r_full = r_inner_full,
      r_notched = r_inner_notched,
      notch_limit = arc_r_inner,
      side = "inner"
    )

    if (nrow(notch_outer) > 0) notch_outer$group <- 1
    if (nrow(notch_inner) > 0) notch_inner$group <- 2

    df_notch_cutout <- bind_rows(notch_outer, notch_inner)

    df_notch_guides <- data.frame(
      x = c(
        r_outer_notched[i1] * cos(notch_lower_plot),
        arc_r_outer * cos(notch_lower_plot),
        r_outer_notched[i2] * cos(notch_upper_plot),
        arc_r_outer * cos(notch_upper_plot),
        r_inner_notched[i1] * cos(notch_lower_plot),
        arc_r_inner * cos(notch_lower_plot),
        r_inner_notched[i2] * cos(notch_upper_plot),
        arc_r_inner * cos(notch_upper_plot)
      ),
      y = c(
        r_outer_notched[i1] * sin(notch_lower_plot),
        arc_r_outer * sin(notch_lower_plot),
        r_outer_notched[i2] * sin(notch_upper_plot),
        arc_r_outer * sin(notch_upper_plot),
        r_inner_notched[i1] * sin(notch_lower_plot),
        arc_r_inner * sin(notch_lower_plot),
        r_inner_notched[i2] * sin(notch_upper_plot),
        arc_r_inner * sin(notch_upper_plot)
      ),
      group = c(1, 1, 2, 2, 3, 3, 4, 4)
    )
  }

  # CI interval is centered within the black box ribbon.
  box_mid <- (median_r_inner + median_r_outer) / 2
  ci_r_inner <- box_mid - ci_thickness / 2
  ci_r_outer <- box_mid + ci_thickness / 2

  # Return every geometry component needed for drawing.
  list(
    side = side,
    theta_plot = theta_plot,
    r_outer_full = r_outer_full,
    r_inner_full = r_inner_full,
    r_outer_notched = r_outer_notched,
    r_inner_notched = r_inner_notched,
    r_outer_out = r_outer_out,
    r_inner_out = r_inner_out,
    df_poly_in = df_poly_in,
    df_poly_out = df_poly_out,
    df_arc_Q1 = arc_ribbon_df(Q1_plot, Median_plot, arc_r_inner, arc_r_outer),
    df_arc_Q3 = arc_ribbon_df(Median_plot, Q3_plot, arc_r_inner, arc_r_outer),
    df_ci = arc_ribbon_df(notch_lower_plot, notch_upper_plot, ci_r_inner, ci_r_outer),
    df_median = arc_ribbon_df(
      Median_plot - 0.025,
      Median_plot + 0.025,
      median_r_inner,
      median_r_outer
    ),
    df_notch_cutout = df_notch_cutout,
    df_notch_guides = df_notch_guides,
    inlier_angles_plot = inlier_angles_plot,
    outlier_angles_plot = outlier_angles_plot,
    statistics = c(
      Q1 = Q1_plot,
      Median = Median_plot,
      Q3 = Q3_plot
    ),
    outliers = outlier_angles_plot
  )
}

