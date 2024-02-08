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
  # NAs/non NAs distinguish valid locations
  origin <- terra::classify(origin, cbind(FALSE, NA))
  # terra::plot(origin)
  origin
}

#' @title Depth error envelopes 

# Depth errors
calc_depth_envelope <- function(.particles, .obs = NULL, .t, .dlist) {
  if (.t == 1L) {
    check_dlist(.dlist = .dlist, 
                .algorithm = "ewindow")
  }
  cell_now <- scale_box <- scale_tail <- NULL
  .particles[, deep := terra::extract(.dlist$algorithm$ewindow$deep, cell_now)]
  .particles[, shallow := terra::extract(.dlist$algorithm$ewindow$shallow, cell_now)]
  .particles[, shallower := terra::extract(.dlist$algorithm$ewindow$shallower, cell_now)]
  .particles[, scale_box := deep - shallow]
  .particles[, scale_tail := shallow - shallower]
 .particles
}

#' @title Depth likelihood term

pf_lik_dc_2 <- function(.particles, .obs, .t, .dlist, .drop) {
  # Checks (one off)
  if (.t == 1L) {
    check_names(.obs, "depth")
    patter:::check_dlist(.dlist = .dlist, 
                         .spatial = c("bathy", "bset"),
                         .algorithm = "ewindow")
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
  .particles <- calc_depth_envelope(.particles = .particles, 
                                    .obs = NULL,
                                    .t = .t, 
                                    .dlist = .dlist)
  # Define depth likelihoods
  # * Likelihood is zero in cells where the depth observation 
  # * ... is not contained within the error window for that location 
  lik_dc <- rep(0L, nrow(.particles))
  pos <- which((depth >= .particles$shallow) & (depth <= .particles$deep))
  if (length(pos) > 0L) {
    lik_dc[pos] <- 0.99 * .particles$scale_box[pos]
  }
  pos <- which((depth >= .particles$shallower) & (depth < .particles$shallow))
  if (length(pos) > 0L) {
    lik_dc[pos] <- 0.01 * .particles$scale_tail[pos]
  }
  # Update likelihoods 
  lik <- NULL
  .particles[, lik := lik * lik_dc, ]
  if (.drop) {
    .particles <- .particles[lik > 0, ]
  }
  .particles
}

#' @title Plot the depth-error model

add_depth_error_model <- function(bathy) {
  stopifnot(inherits(bathy, "numeric"))
  stopifnot(length(bathy) == 1L)
  # Define window of possible depth observations given a bathymetric depth 
  deep           <- edeep(.bathy = bathy)
  shallow        <- eshallow(.bathy = bathy)
  shallower_howe <- eshallower(.eshallow = shallow, .emove = emovehowe)
  shallower_digi <- eshallower(.eshallow = shallow, .emove = emovedigi)
  # Scale depth limits
  if (FALSE) {
    deep           <- deep - bathy
    shallow        <- shallow - bathy 
    shallower_howe <- shallower_howe - bathy
    shallower_digi <- shallower_digi - bathy
  }
  # Add lines to a plot
  # * TO DO
  # * Update scaling
  lines(c(shallower_digi, shallower_digi), c(0, 0.01), col = "red")
  lines(c(shallower_howe, shallower_howe), c(0, 0.01))
  lines(c(shallower_digi, shallower_howe), c(0.01, 0.01), col = "red")
  lines(c(shallower_howe, shallow), c(0.01, 0.01))
  lines(c(shallow, shallow), c(0.01, 1))
  lines(c(shallow, deep), c(1, 1))
  lines(c(deep, deep), c(1, 0))
}

#' @title Depth likelihood model for coarse grid 
#' (see R/supporting/trial-coarse-fine-2.R)
#' TO DO
# * Review `cell_now` use is correct

pf_lik_dc_coarse <- function(.particles, .obs, .t, .dlist, drop) {
  
  if (.t == 1L) {
    check_dlist(.dlist = .dlist, 
                .spatial = c("bathy", "original", "v"))
  }
  
  #### Identify unique grid cells 
  pu <- 
    .particles |> 
    select("cell_now", "x_now", "y_now") |>
    distinct(.data$cell_now, .keep_all = TRUE) |>
    as.data.table()
  
  #### Define 'buffers'
  buffer <- 
    .dlist$spatial$v[pu$cell_now] |> 
    sf::st_as_sf()
  
  #### Define likelihoods
  # Use data.table method
  # For each coarse grid cell, extract all 'internal' cells
  # & ... filter invalid ones & identify thus identify valid coarse cells
  ff <- function(df, ...) {
    # Evaluate likelihood of (internal) cells
    pnow <- 
      df |> 
      select(cell_now = "cell", x_now = "x", y_now = "y", bathy = "value") |> 
      mutate(lik = 1L) |>
      as.data.table() |>
      pf_lik_dc_2(.obs = .obs, .t = .t, .dlist = .dlist, drop = .drop) 
    # Define whether or not there are any valid (internal) cells
    nrow(pnow) > 0
  }
  # Identify valid proposals
  valid <- NULL
  pu[, valid := exactextractr::exact_extract(.dlist$spatial$original, 
                                             buffer, 
                                             include_cell = TRUE, 
                                             include_xy = TRUE,
                                             summarize_df = TRUE, 
                                             fun = ff, 
                                             progress = TRUE) |> unlist()]
  
  #### Return valid proposals
  cell_now <- NULL
  .particles[cell_now %in% pu$cell_now[pu$valid], ]
}
