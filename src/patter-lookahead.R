#' @title Normalise a numeric vector

normalise <- function(.x) {
  .x / sum(.x)
}

#' @title Planar coordinates' centroid 

geocentroid <- function(.xy) {
  matrix(apply(.xy, 2, mean), nrow = 1L, ncol = 2L)
}

#' @title Lookahead model

pf_lik_ac_lookahead <- function(.particles, .obs, .t, .dlist) {

  #### Check user inputs
  check_dlist(.dlist = .dlist,
              .algorithm = "dlen",
              .par = c("lonlat", "shape", "scale", "mobility", "gamma"))

  #### Define time gap before detection
  # Define time gap
  pos_detections <- .dlist$algorithm$pos_detections
  pos_detection  <- (pos_detections[pos_detections > .t])[1L]
  timegap <- pos_detection - .t
  # For long time gaps, we skip this step
  # * Euclidean distances are not very meaningful
  # * We do not have to worry about convergence until we get nearer to a receiver
  if (timegap > 6 * 30L) {
    return(.particles) 
  }
  
  #### Define anchor point (average receiver coordinate)
  # We treat this as a single point of 'highest likelihood'
  axy <-
    .dlist$data$moorings |>
    filter(.data$receiver_id %in% .obs$receiver_id_next[.t][[1]]) |>
    select("receiver_x", "receiver_y") |>
    as.data.table() |>
    as.matrix() |>
    geocentroid()

  #### Calculate distance from particles to anchor point
  dist <- terra::distance(.particles[, list(x_now, y_now)] |>
                            as.matrix(),
                          axy,
                          lonlat = .dlist$pars$lonlat)

  #### Define parameters for density calculations 
  # Extract parameters for one time step
  sh <- .dlist$pars$shape
  sc <- .dlist$pars$scale
  ga <- .dlist$pars$gamma
  mb <- .dlist$pars$mobility
  # Update parameters in line with time gap
  # (and accounting for size of detection range relative to mobility)
  timegap <- timegap + ga / mb 
  sh <- sh * timegap
  mb <- mb * timegap
  
  #### Translate distances into probability densities
  dens <- .dlist$algorithm$dlen(.x = dist, .shape = sh, .scale = sc, .mobility = mb)

  #### Update 'likelihood' term (used for resampling)
  # * This eliminates particles that are too far away from from the next receiver
  # * It also upweights particles that are nearer the next receiver, facilitating convergence
  lik <- NULL
  .particles[, lik := normalise(lik * dens)]
  .particles
}