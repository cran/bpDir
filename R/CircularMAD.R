#' @export
#' @importFrom graphics hist
# Circular Median Absolute Deviation (CMAD)
CircularMAD <- function(X, h = FALSE) {

  # Compute circular median
  SFMAngle <- median.circular(circular(X, modulo = "2pi"))

  # Convert points and median to Cartesian coordinates
  X_cart <- pol2cart_circle(X)
  SFM <- pol2cart_circle(SFMAngle)

  # Arc distances between median and points
  arcdist <- acos(SFM %*% t(X_cart))

  # Optional histogram visualization
  if (h == TRUE) {
    hist(arcdist, breaks = seq(from = 0, to = pi, by = pi/8))
    points(c(0) ~ median(arcdist), col = "blue", pch = 19)
    points(c(0) ~ mean(arcdist), col = "red", pch = 19)
  }

  # Return CMAD
  median(arcdist)
}
