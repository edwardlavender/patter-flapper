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
ebathy <- function(.bathy, .a = 0.5, .b = 0.013) {
  sqrt(.a + (.b * .bathy)^2)
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

# Simulates time series of behavioural states
simulate_behaviour <- function(timeline) {
  
  # Simulate initial state (resting/active)
  # > Individuals spend approx 50 % time resting on average
  beh_state   <- sample(c(1L, 2L), size = 1)
  
  # Simulate subsequent states
  beh_ts      <- list()
  chunk <- 1L
  span  <- 1
  while (span < length(timeline)) {
    
    # Simulate the number of time steps for which the current state is maintained
    # * These parameters are based on explore-data-mefs-behaviour.R
    if (beh_state == 1L) {
      z <- rgamma(n = 1L, shape = 0.730, rate = 0.163)
    } else {
      z <- rgamma(n = 1L, shape = 0.960, rate = 0.168)
    }
    z <- plyr::round_any(z, 2)
    
    # Record duration
    beh_ts[[chunk]] <- data.table(state = beh_state, duration = z)
    
    # Update loop
    chunk     <- chunk + 1L
    beh_state <- ifelse(beh_state == 1L, yes = 2L, no = 1L)
    span      <- sum(rbindlist(beh_ts)$duration)
    
  }
  
  # Define time series of states
  beh_ts <- rbindlist(beh_ts)
  beh_ts <- beh_ts[rep(1:.N, duration), ]
  beh_ts <- beh_ts[1:length(timeline), ]
  stopifnot(nrow(beh_ts) == length(timeline))
  stopifnot(all(!is.na(beh_ts$state)))
  
  # Examine properties of simulated dataset
  # > Approximately 50 % of the time spent resting 
  if (FALSE) {
    table(beh_ts$state)
    duration_state(beh_ts, fct = NULL, state = 0L) |> pull(duration) |> hist()
    duration_state(beh_ts, fct = NULL, state = 1L) |> pull(duration) |> hist()
  }
  
  # Return sequence of states
  beh_ts$state
  
}

# Truncated gamma distributions for step length 
rtruncgamma <- function(.n, .shape, .scale, .mobility) {
  truncdist::rtrunc(.n, "gamma", a = 0, b = .mobility,
                    shape = .shape, scale = .scale)
}

dtruncgamma <- function(.x, .shape, .scale, .mobility) {
  truncdist::dtrunc(.x, "gamma", a = 0, b = .mobility,
                    shape = .shape, scale = .scale)
}

