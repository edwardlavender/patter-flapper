#' @title Depth-error model parameters

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
  # terra::plot(origin)
  origin
}

#' @title Depth error envelopes 

# Depth errors
calc_depth_envelope <- function(.particles, .obs = NULL, .t, .dlist) {
  if (.t == 1L) {
    check_names(.obs, c("bathy", "digi"))
    check_dlist(.dlist = .dlist, 
                .algorithm = "ewindow")
  }
  cell_now <- NULL
  particles[, deep := terra::extract(.dlist$algorithm$ewindow$deep, cell_now)]
  particles[, shallow := terra::extract(.dlist$algorithm$ewindow$deep, cell_now)]
  particles[, shallower := terra::extract(.dlist$algorithm$ewindow$deep, cell_now)]
 .particles
}

#' @title Depth likelihood term

pf_lik_dc_2 <- function(.particles, .obs, .t, .dlist) {
  # Checks (one off)
  if (.t == 1L) {
    patter:::check_dlist(.dlist = .dlist, 
                         .spatial = c("bathy", "bset"),
                         .algorithm = c("ewindow"))
  }
  # Handle NA observations
  depth <- .obs$depth[.t]
  if (is.na(depth)) {
    return(.particles)
  }
  # Define bathymetric depth & data source (digi versus howe)
  cell_now <- NULL
  if (!rlang::has_name(.particles, "bathy")) {
    .particles[, bathy := terra::extract(.dlist$spatial$bathy, cell_now)]
  }
  .particles[, digi := terra::extract(.dlist$spatial$bset, cell_now)]
  # Define depth envelope
  .particles <- calc_depth_envelope(.particles = .particles)
  # Define depth likelihoods
  # * Likelihood is zero in cells where the depth observation 
  # * ... is not contained within the error window for that location 
  lik_dc <- rep(0L, nrow(.particles))
  pos <- which((depth >= .dlist$algorithm$ewindow$shallow) & (depth <= .dlist$algorithm$ewindow$deep))
  lik_dc[pos] <- 0.99
  pos <- which((depth >= .dlist$algorithm$ewindow$shallower) & (depth <= .dlist$algorithm$ewindow$deep))
  lik_dc[pos] <- 0.01
  # Update likelihoods 
  lik <- NULL
  .particles[, lik := lik * lik_dc, ]
  .particles[lik > 0, ]
}

#' @title Plot the depth-error model

add_depth_error_model <- function(bathy) {
  stopifnot(inherits(bathy, "numeric"))
  stopifnot(length(bathy) == 1L)
  # Define window of possible depth observations given a bathymetric depth 
  deep           <- edeep(.bathy = bathy)
  shallow        <- eshallow(.bathy = bathy)
  shallower_howe <- eshallower(.eshallow = bathy, .emove = 0L)
  shallower_digi <- eshallower(.eshallow = bathy, .emove = 1L)
  # Scale depth limits
  if (TRUE) {
    deep           <- deep - bathy
    shallow        <- shallow - bathy 
    shallower_howe <- shallower_howe - bathy
    shallower_digi <- shallower_digi - bathy
  }
  # Add lines to a plot
  lines(c(shallower_digi, shallower_digi), c(0, 0.01), col = "red")
  lines(c(shallower_howe, shallower_howe), c(0, 0.01))
  lines(c(shallower_digi, shallower_howe), c(0.01, 0.01), col = "red")
  lines(c(shallower_howe, shallow), c(0.01, 0.01))
  lines(c(shallow, shallow), c(0.01, 1))
  lines(c(shallow, deep), c(1, 1))
  lines(c(deep, deep), c(1, 0))
}
