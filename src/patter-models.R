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

# Movement error (into shallower waters)
# * For some receivers in the Howe data, the individual is ~10-20 m shallower than previously assumed.
# * For the Digimap data, the individual can be up to 40 m shallower. 
emovehowe <- 20
emovedigi <- 40
emove <- function(.digi) {
  e <- data.frame(digi = c(0L, 1L), 
                  error = c(emovehowe, emovedigi))
  e$error[match(.digi, e$digi)]
}

# Shallow and deep limits for observed depth in given location
# * .bathy may be a SpatRaster or a vector of values
# * other terms are single numbers
eshallow <- function(.bathy, 
                     .ebathy = ebathy(.bathy), 
                     .etag = etag, 
                     .etide = etide) {
  .bathy - .ebathy - .etag - .etide - 5
}

eshallower <- function(.eshallow,
                       .emove) {
  # .bathy - .ebathy - .etag - .etide - .emove
  .eshallow - .emove
}

edeep <- function(.bathy, 
                  .ebathy = ebathy(.bathy), 
                  .etag = etag, 
                  .etide = etide) {
  .bathy + .ebathy + .etag + .etide
}

#' @title Origin spatRaster for *DC models (DCPF, ACDCPF)

# Define the min/max possible depth of the individual in each location
# * This is independent of the depth observation
# * We then assess the likelihood of the depth observation in each cell
# * ... based on whether or not it overlaps with the required window
dc_ewindow <- function(.bathy, .bset) {
  # Define spatially explicit error terms 
  # * (80 s + 12 s + 12 s)
  ebathyspat <- ebathy(.bathy)
  emovespat  <- terra::classify(.bset, cbind(0, emovehowe))
  emovespat  <- terra::classify(emovespat, cbind(1, emovedigi))
  # Define possible depth range for individual in each location 
  # * (65 s + 25 s + 45 s)
  shallow   <- eshallow(.bathy = .bathy, 
                        .ebathy = ebathyspat)
  shallower <- eshallower(.eshallow = shallow, 
                          .emove = emovespat)
  deep      <- edeep(.bathy = .bathy, 
                     .ebathy = ebathyspat)
  # Return list
  list(shallow = shallow, shallower = shallower, deep = deep)
}

# Define possible locations of the individual on a SpatRaster
# * .ewindow is the depth window list
# * .depth is the observed depth 
dc_origin <- function(.ewindow, .depth) {
  # ~53 s
  # * Likelihood = 1 in cells where the depth-error window contains the observed depth
  # * Likelihood = 0 in cells where the depth-error window doesn't contain the observed depth 
  origin <- (.depth >= .ewindow$shallower) & (.depth <= .ewindow$deep)
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