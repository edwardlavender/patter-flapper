#' Acoustic observations

ddetlogistic <- function(.x, .alpha, .beta, .gamma) {
  ifelse(.x <= .gamma, plogis(.alpha + .beta * .x), 0)
}

#' Archival observations 

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

#' Movement models

rtruncgamma <- function(.n, .shape, .scale, .mobility) {
  truncdist::rtrunc(.n, "gamma", a = 0, b = .mobility,
                    shape = .shape, scale = .scale)
}

dtruncgamma <- function(.x, .shape, .scale, .mobility) {
  truncdist::dtrunc(.x, "gamma", a = 0, b = .mobility,
                    shape = .shape, scale = .scale)
}