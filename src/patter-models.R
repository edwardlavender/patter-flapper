#' @title Acoustic observations

ddetlogistic <- function(.x, .alpha, .beta, .gamma) {
  ifelse(.x <= .gamma, plogis(.alpha + .beta * .x), 0)
}

#' @title Archival observations

# Tag error
etag   <- 4.77

# Tidal range 
etide  <- 2.50

# Bathymetric error (depth-dependent)
# * .bathy may be a SpatRaster or a numeric vector
ebathy <- function(.bathy) {
  sqrt(0.5002 + (0.013 * .bathy)^2)
}

# Total deep error
# (See also process-data-spatial-raster.R for full justification)
# ebathy + etag + etide + extra
ebathy(350) + etag + etide
edeep <- 20

# Total shallow area 
# ebathy + etag + etide + edemersal (5 m) + extra
ebathy(350) + etag + etide + 5
eshallow <- 20

#' @title Origin spatRaster for *DC models (DCPF, ACDCPF)

dc_origin <- function(.map, .depth) {
  # Calculate depth window
  deep    <- .map - edeep
  shallow <- .map + eshallow
  # Evaluate 'likelihood' (~53 s)
  # * Likelihood = 1 in cells where the depth-error window contains the observed depth
  # * Likelihood = 0 in cells where the depth-error window doesn't contain the observed depth 
  origin <- (.depth >= shallow) & (.depth <= deep)
  # NAs/non NAs distinguish valid locations
  origin <- terra::classify(origin, cbind(FALSE, NA))
  # terra::plot(origin)
  origin
}

#' @title Movement models

rtruncgamma <- function(.n, .shape, .scale, .mobility) {
  truncdist::rtrunc(.n, "gamma", a = 0, b = .mobility,
                    shape = .shape, scale = .scale)
}

dtruncgamma <- function(.x, .shape, .scale, .mobility) {
  truncdist::dtrunc(.x, "gamma", a = 0, b = .mobility,
                    shape = .shape, scale = .scale)
}