# circular_violin_plot()
#
# This function draws circular violin plots with options:
#   - full single, outer single, inner single, or two-halves inner outer
#   - notches
#   - bootstrap confidence interval
#   - inlier and outlier density regions
#   - beans for observations or outlier points
#
# The notch is intentionally drawn before the black box so the box masks
# the end of the notch lines and prevents boundary artifacts.
#
# ------------------------------------------------------------
#' @export
#' @importFrom ggplot2 ggplot geom_path aes geom_polygon geom_point coord_equal theme_void labs scale_linewidth_identity

CircularViolinPlot <- function(
    data,
    data2 = NULL,
    type = c("full", "outer", "inner", "two_halves"),
    notch = FALSE,
    beans = c("none", "all", "inliers", "outliers", "outlier_points"),
    outlier_point_color = "red",
    outlier_point_size = 2,
    outlier_point_alpha = 0.8,
    show_ci_interval = TRUE,
    show_notch_cutout = TRUE,
    show_notch_guides = FALSE,
    r0 = 1,
    density_scale = 0.3,
    vonMises_grid = NULL,
    bw = NULL,
    bw2 = NULL,  # Optional bandwidth for data2 when type = "two_halves"
    kernel = "vonmises",
    fill = "skyblue",
    fill2 = "darkorange",
    alpha = 0.75,
    alpha2 = 0.75,
    fill_outlier = NULL,
    fill2_outlier = NULL,
    alpha_outlier = 0.2,
    alpha2_outlier = 0.2,
    arc_color = "black",
    arc_color2 = "black",
    ci_color = "white",
    ci_color2 = "white",
    arc_thickness = 0.06,
    ci_thickness = 0.012,
    box_width_factor = 1,
    density_threshold = 0,
    cutoff = NULL,
    cutoff2 = NULL,  # Optional cutoff for data2 when type = "two_halves"
    #zero = 0,
    #rotation = "counter",
    n_boot = 1000,
    ci_level = 0.95,
    notch_factor = 0.9,
    notch_guide_linetype = "solid",
    notch_guide_linewidth = 0.4,
    bean_length = 0.035,
    bean_min_size = 0.15,
    bean_max_size = 0.6,
    circle_color = "grey70",
    circle_linewidth = 0.7,
    center_point = TRUE,
    center_point_color = "black",
    center_point_size = 1.8,
    title = "",
    xlim = c(-1.5, 1.5),
    ylim = c(-1.5, 1.5)
) {

  # ----------------------------------------------------------
  # Load required packages.
  # circular: circular statistics and circular density estimation.
  # ggplot2: plotting.
  # dplyr: bind_rows() and data-frame utilities.
  # ----------------------------------------------------------
  #library(circular)
  #library(ggplot2)
  #library(dplyr)

  type <- match.arg(type) #default to 'full' if not specified by user
  beans <- match.arg(beans)

  if(!is.circular(data)){stop("argument data must be entered as a vector of class circular")}

  #Record properties of input vector, to use later in the output
  typeOfData <- circular::circularp(data)$type
  unitsOfData <- circular::circularp(data)$units
  templateOfData <- circular::circularp(data)$template
  moduloOfData <- circular::circularp(data)$modulo
  zeroOfData <- circular::circularp(data)$zero
  rotationOfData <- circular::circularp(data)$rotation


  data  <- conversion.circular(data, units = "radians", modulo="2pi", zero=0, rotation="counter")

  if (!is.null(data2)) {

    if(!is.circular(data2)){stop("argument data2 must be entered as a vector of class circular")}

    #Record properties of input vector, to use later in the output
    typeOfData2 <- circular::circularp(data)$type
    unitsOfData2 <- circular::circularp(data)$units
    moduloOfData2 <- circular::circularp(data)$modulo
    templateOfData2 <- circular::circularp(data2)$template
    zeroOfData2 <- circular::circularp(data)$zero
    rotationOfData2 <- circular::circularp(data)$rotation

    data2  <- conversion.circular(data2, units = "radians", modulo="2pi", zero=0, rotation="counter")

  }






  if (is.null(bw)) {

    if (is.null(vonMises_grid)) {
      gvm <- generate_vonmises_grid()
    } else {
      gvm <- vonMises_grid
    }

    kclsm <- EstimateKappa(data, "CLMS", gvm)

    bw <- kclsm$kappa

    cat("bw = \n")
    print(bw)

  }

  if (is.null(cutoff)) {

    if (is.null(vonMises_grid)) {
      gvmForCutoff <- generate_vonmises_grid()
    } else {
      gvmForCutoff <- vonMises_grid
    }

    kclsmForCutoff <- EstimateKappa(data, "CLMS", gvmForCutoff)

    cutoff <-qvonmises(.995, mu=circular(0), kappa=kclsmForCutoff$kappa)

    cat("cutoff = \n")
    print(conversion.circular(cutoff, type = typeOfData, units = unitsOfData, template = templateOfData, modulo= moduloOfData, zero=zeroOfData, rotation=rotationOfData))

  }

  # Validate user choices for `type` and `beans`.
  #type <- match.arg('type')
  #beans <- match.arg('beans')

  # two_halves needs two datasets: data for the outer half and data2 for the inner half.
  if (type == "two_halves" && is.null(data2)) {
    stop("`data2` is required when type = 'two_halves'.")
  }

  # Set defaults that depend on other arguments.
  # cutoff = pi gives a half-circle fence around the median.
  # If outlier fill colors are not supplied, reuse the corresponding main fill.
  if (is.null(cutoff)) cutoff <- pi
  # By default, data2 uses the same cutoff as data unless cutoff2 is supplied.
  if (is.null(cutoff2)) cutoff2 <- cutoff
  if (is.null(fill_outlier)) fill_outlier <- fill
  if (is.null(fill2_outlier)) fill2_outlier <- fill2

  # ----------------------------------------------------------
  # Build geometry according to requested plot type.
  # ----------------------------------------------------------
  if (type == "full") {
    geoms <- list(build_geom(data,n_boot,ci_level,kernel,density_scale,density_threshold,notch,notch_factor,r0,arc_thickness,box_width_factor,ci_thickness, side="both", bw_use_arg = bw, cutoff_use_arg = cutoff))
    fills <- fill
    alphas <- alpha
    fill_outliers <- fill_outlier
    alpha_outliers <- alpha_outlier
    arcs <- arc_color
    cis <- ci_color
  }

  if (type == "outer") {
    geoms <- list(build_geom(data,n_boot,ci_level,kernel,density_scale,density_threshold,notch,notch_factor,r0,arc_thickness,box_width_factor,ci_thickness, side = "outer", bw_use_arg = bw, cutoff_use_arg = cutoff))
    fills <- fill
    alphas <- alpha
    fill_outliers <- fill_outlier
    alpha_outliers <- alpha_outlier
    arcs <- arc_color
    cis <- ci_color
  }

  if (type == "inner") {
    geoms <- list(build_geom(data,n_boot,ci_level,kernel,density_scale,density_threshold,notch,notch_factor,r0,arc_thickness,box_width_factor,ci_thickness, side="inner", bw_use_arg = bw, cutoff_use_arg = cutoff))
    fills <- fill
    alphas <- alpha
    fill_outliers <- fill_outlier
    alpha_outliers <- alpha_outlier
    arcs <- arc_color
    cis <- ci_color
  }

  if (type == "two_halves") {

    if (is.null(bw2)) {


      if (is.null(vonMises_grid)) {
        gvm2 <- generate_vonmises_grid()
      } else {
        gvm2 <- vonMises_grid
      }

      kclsm2 <- EstimateKappa(data2, "CLMS", gvm2)

      bw2 <- kclsm2$kappa

      cat("bw2 = \n")
      print(bw2)

    }

    if (is.null(cutoff)) {

      if (is.null(vonMises_grid)) {
        gvmForCutoff2 <- generate_vonmises_grid()
      } else {
        gvmForCutoff2 <- vonMises_grid
      }

      kclsmForCutoff2 <- EstimateKappa(data2, "CLMS", gvmForCutoff2)

      cutoff2 <-qvonmises(.995, mu=circular(0), kappa=kclsmForCutoff2$kappa)

      cat("cutoff2 = \n")
      print(conversion.circular(cutoff2, type = typeOfData2, units = unitsOfData2, template = templateOfData2, modulo= moduloOfData2, zero=zeroOfData2, rotation=rotationOfData2))

    }

    geoms <- list(
      build_geom(data, n_boot,ci_level,kernel,density_scale,density_threshold,notch,notch_factor,r0,arc_thickness,box_width_factor,ci_thickness, side="outer", bw_use_arg = bw, cutoff_use_arg = cutoff),
      build_geom(data2,n_boot,ci_level,kernel,density_scale,density_threshold,notch,notch_factor,r0,arc_thickness,box_width_factor,ci_thickness, side="inner", bw_use_arg = bw2, cutoff_use_arg = cutoff2)
    )
    fills <- c(fill, fill2)
    alphas <- c(alpha, alpha2)
    fill_outliers <- c(fill_outlier, fill2_outlier)
    alpha_outliers <- c(alpha_outlier, alpha2_outlier)
    arcs <- c(arc_color, arc_color2)
    cis <- c(ci_color, ci_color2)
  }

  #make all plotting geoms available in output
  #return(list(
  #  geoms = geoms
  #))

  # Reference unit circle coordinates.
  theta_seq <- seq(0, 2*pi, length.out = 500)
  df_circle <- data.frame(
    x = cos(theta_seq),
    y = sin(theta_seq)
  )

  # Start with an empty ggplot.
  p <- ggplot()

  # Draw the grey reference circle.
  p <- p +
    geom_path(
      data = df_circle,
      aes(x, y),
      color = circle_color,
      linewidth = circle_linewidth
    )

  # Add all geometry objects to the plot.
  # two_halves has two geometry objects; other types have one.
  for (i in seq_along(geoms)) {
    p <- add_geom(
      p,
      geoms[[i]],
      fill = fills[i],
      alpha = alphas[i],
      fill_outlier = fill_outliers[i],
      alpha_outlier = alpha_outliers[i],
      arc_color = arcs[i],
      ci_color = cis[i],
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
    )
  }

  # Use identity linewidth scale because bean linewidths are already computed.
  if (beans %in% c("all", "inliers", "outliers")) {
    p <- p + scale_linewidth_identity()
  }


  # Draw CI interval after the grey reference circle so it remains visible
  # in full, outer, inner, and two_halves modes.
  if (show_ci_interval) {
    for (i in seq_along(geoms)) {
      if (nrow(geoms[[i]]$df_ci) > 0) {
        p <- p +
          geom_polygon(
            data = geoms[[i]]$df_ci,
            aes(x, y),
            fill = cis[i],
            color = cis[i],
            linewidth = 0.15
          )
      }
    }
  }

  # Optional center mark.
  if (center_point) {
    p <- p +
      geom_point(
        aes(x = 0, y = 0),
        shape = 4,
        color = center_point_color,
        size = center_point_size,
        stroke = 1
      )
  }

  # Final plot settings.
  p <- p +
    coord_equal(xlim = xlim, ylim = ylim) +
    theme_void() +
    labs(title = title)

  if (type == "two_halves") {
    return(list(
      plot = p,
      statistics = list(
        data = conversion.circular(geoms[[1]]$statistics, type = typeOfData, units = unitsOfData, template = templateOfData, modulo= moduloOfData, zero=zeroOfData, rotation=rotationOfData),
        data2 = conversion.circular(geoms[[2]]$statistics, type = typeOfData2, units = unitsOfData2, template = templateOfData2, modulo= moduloOfData2, zero=zeroOfData2, rotation=rotationOfData2)
      ),
      outliers = list(
        data = conversion.circular(geoms[[1]]$outliers, type = typeOfData, units = unitsOfData, template = templateOfData, modulo= moduloOfData, zero=zeroOfData, rotation=rotationOfData),
        data2 = conversion.circular(geoms[[2]]$outliers, type = typeOfData2, units = unitsOfData2, template = templateOfData2, modulo= moduloOfData2, zero=zeroOfData2, rotation=rotationOfData2)
      )
    ))
  }

  return(list(
    plot = p,
    statistics = conversion.circular(geoms[[1]]$statistics, type = typeOfData, units = unitsOfData,template = templateOfData, modulo= moduloOfData, zero=zeroOfData, rotation=rotationOfData),
    outliers = conversion.circular(geoms[[1]]$outliers,  type = typeOfData, units = unitsOfData, template = templateOfData, modulo= moduloOfData, zero=zeroOfData, rotation=rotationOfData)
  ))

}

