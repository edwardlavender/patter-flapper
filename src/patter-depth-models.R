#' @title Depth error function
# https://github.com/edwardlavender/flapper_appl/blob/master/R/define_global_param.R

calc_depth_error <- readRDS(dv::here_data("input", "pars.rds"))$patter$calc_depth_error

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
  e <- data.frame(digi = c(0, 1), 
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
  (.bathy - .ebathy) - (.etag - .etide) - 5
}
eshallower <- function(.eshallow,
                       .emove) {
  # (.bathy - .ebathy) - (.etag - .etide) - (.emove)
  .eshallow - .emove
}
edeep <- function(.bathy, 
                  .ebathy = ebathy(.bathy), 
                  .etag = etag, 
                  .etide = etide) {
  (.bathy + .ebathy) + (.etag + .etide)
}

#' @title Origin spatRaster for *DC models (DCPF, ACDCPF)

dc_origin <- function(.bathy, .bset, .depth) {
  
  # Define spatially explicit error terms
  ebathyspat <- ebathy(.bathy)
  emovespat  <- terra::classify(.bset, cbind(0, emovehowe))
  emovespat  <- terra::classify(emovespat, cbind(1, emovedigi))
  
  # Define possible depth range for individual in each location
  shallow   <- eshallow(.bathy = .bathy, 
                        .ebathy = ebathyspat)
  shallower <- eshallower(.eshallow = shallow, 
                          .emove = emovespat)
  deep      <- edeep(.bathy = .bathy, 
                     .ebathy = ebathyspat)

  # Define regions on bathy 
  # * Regions beyond the shallow and deep limits are NA
  # * The likelihood of the depth observation in these cells is 0
  .bathy <- terra::mask(.bathy, .bathy < shallower, maskvalues = TRUE)
  .bathy <- terra::mask(.bathy, .bathy > deep, maskvalues = TRUE)
  # terra::plot(.bathy)
  .bathy
}

#' @title Depth error envelopes 

# Depth errors
calc_depth_envelope <- function(.particles) {
  # Requirements:
  # * bathy column
  # * digi column
  bathy <- digi <- NULL
  .particles[, ebathyterm := edeep(.bathy = bathy)]
  .particles[, deep := edeep(.bathy = bathy, .ebathy = ebathyterm)]
  .particles[, shallow := eshallow(.bathy = bathy, .bathy = ebathyterm)]
  .particles[, shallower := eshallower(.eshallow = shallow, .emove = emove(digi))]
 .particles
}

#' @title Depth likelihood term

pf_lik_dc_2 <- function(.particles, .obs, .t, .dlist) {
  # Checks (one off)
  if (.t == 1L) {
    patter:::check_dlist(.dlist = .dlist, 
                         .spatial = "bathy", 
                         .algorithm = c("bset", "calc_depth_error"))
  }
  # Define bathymetric depth & data source (digi versus howe)
  cell_now <- NULL
  if (!rlang::has_name(.particles, "bathy")) {
    .particles[, bathy := terra::extract(.dlist$spatial$bathy, cell_now)]
  }
  .particles[, digi := terra::extract(.dlist$algorithm$bset, cell_now)]
  # Define depth envelope
  .particles <- calc_depth_envelope(.particles = .particles)
  # Define depth likelihoods
  # * Likelihood is zero in cells whether bathymetric depth does not fall between required limits
  lik_dc <- rep(0L, nrow(.particles))
  pos <- which(.particles$bathy <= .particles$deep & 
                 .particles$bathy  >= .particles$shallow)
  lik_dc[pos] <- 0.99
  pos <- which(.particles$bathy < .particles$shallow & 
                 .particles$bathy >= .particles$shallower)
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
