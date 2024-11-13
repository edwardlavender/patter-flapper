#' Define a 2 minute timeline given the mmyy
patter_timeline <- function(mmyy) {
  period <- mmyyrng(mmyy)
  seq(period[1],  period[2], by = "2 mins")
}

#' Define the model_obs and yobs inputs for a forward run of the particle filter
patter_ModelObs <- function(sim, timeline, detections, moorings, archival, xinit_list) {
  
  #### Acoustics datasets
  if (sim$dataset %in% c("ac", "acdc")) {
    # Detection model parameters
    receiver_alpha <- receiver_beta <- receiver_gamma <- NULL
    if (rlang::has_name(moorings, "receiver_alpha")) {
      moorings[, receiver_alpha := NULL]
    }
    if (rlang::has_name(moorings, "receiver_beta")) {
      moorings[, receiver_beta := NULL]
    }
    if (rlang::has_name(moorings, "receiver_gamma")) {
      moorings[, receiver_gamma := NULL]
    }
    moorings[, receiver_alpha := sim$receiver_alpha]
    moorings[, receiver_beta := sim$receiver_beta]
    moorings[, receiver_gamma := sim$receiver_gamma]
    # Assemble acoustic data
    acoustics <- assemble_acoustics(.timeline = timeline,
                                    .detections = detections, 
                                    .moorings = moorings)
    # Assemble container data
    # max(c(ymax(bathy) - ymin(bathy), xmax(bathy) - xmin(bathy)))
    # Note that containers contains a $forward and $backward element that we select later
    acoustics_containers <- assemble_acoustics_containers(.timeline  = timeline, 
                                                          .acoustics = acoustics,
                                                          .mobility  = sim$mobility, 
                                                          .map = qs::qread(here_data("spatial", "bathy-bbox.qs")))
  }
  
  #### Archival datasets
  # This is only required for real-world data
  if (sim$dataset %in% c("dc", "acdc")) {
    if (!all(c("timestamp", "sensor_id", "obs", "depth_sigma", "depth_deep_eps") %in% colnames(archival))) {
      archival <- assemble_archival(.timeline = timeline, 
                                    .archival = 
                                      archival |> 
                                      rename(obs = "depth") |> 
                                      mutate(sensor_id = 1L, 
                                             depth_sigma = sim$depth_sigma, 
                                             depth_deep_eps = sim$depth_deep_eps) |>
                                      select("timestamp", "sensor_id", "obs", "depth_sigma", "depth_deep_eps") |> 
                                      as.data.table()
      )
    } else {
      warn("Using `archival` as inputed, including with existing $depth_sigma and $depth_deep_eps columns.")
    }
  }
  
  #### Capture locations
  # xinit_list is a list with forward & backward elements
  # capture_containers contains a $forward and $backward element that we select later
  capture_containers <- assemble_capture_containers(.timeline  = timeline, 
                                                    .xinit     = xinit_list, 
                                                    .radius    = 500, 
                                                    .mobility  = sim$mobility, 
                                                    .threshold = 139199, 
                                                    .as_ModelObsAcousticContainer = TRUE)
  # (optional) TO DO
  # Improve efficiency of capture containers if ModelObsAcousticContainer is also used
  # We only need capture containers after the last acoustic container (moving forwards in time)

  #### Return list of algorithm-specific ModelObs strings & corresponding datasets
  # Handle containers
  # * Combine acoustic and capture containers for improved efficiency 
  # * (<= 3 ModelObs types are MUCH faster than >= 4 types)
  if (sim$dataset %in% c("ac", "acdc")) {
    containers <- acoustics_containers
    for (direction in c("forward", "backward")) {
      containers[[direction]] <- 
        rbind(containers[[direction]], capture_containers[[direction]]) |> 
        arrange(timestamp, sensor_id) |> 
        as.data.table()
    }
  } else {
    containers <- capture_containers
  }
  # Define list
  if (sim$dataset == "ac") {
    out <- list(ModelObsAcousticLogisTrunc = acoustics, 
                ModelObsAcousticContainer  = containers)
  } else if (sim$dataset == "dc") {
    out <- list(ModelObsDepthNormalTrunc = archival, 
                ModelObsAcousticContainer = containers)
  } else if (sim$dataset == "acdc") {
    out <- list(ModelObsAcousticLogisTrunc = acoustics,
                ModelObsAcousticContainer  = containers, 
                ModelObsDepthNormalTrunc   = archival)
  }
  
  out
  
}

# #' Define the number of particles to use for ACPF/DCPF/ACDCPF runs
# patter_np <- function(sim) {
#   if (sim$dataset == "ac") {
#     np <- 1e4L
#   } else if (sim$dataset == "dc") {
#     np <- 1e5L
#   } else if (sim$dataset == "acdc") {
#     np <- 1e5L
#   }
#   np
# }

#' Define a pf_filter() that runs the algorithms for a given simulation
pf_filter_wrapper <- function(sim, args) {
  
  # browser()
  
  # Define outputs (fwd, bwd)
  if (args$.direction == "forward") {
    routine <- "fwd"
    # pfile   <- "coord-fwd.qs"
    dfile   <- "data-fwd.qs"
  } else if (args$.direction == "backward") {
    routine <- "bwd"
    # pfile   <- "coord-bwd.qs"
    dfile   <- "data-bwd.qs"
  } else {
    stop("args$.direction is missing!")
  }
  
  # Define data & parameters
  # * These are passed down from args, above
  
  # Initialise variables
  success     <- FALSE
  error       <- NA_character_
  n_trial     <- NA_integer_

  # Run filter
  t1      <- Sys.time()
  pout    <- tryCatch(do.call(pf_filter, args, quote = TRUE), error = function(e) e)
  t2      <- Sys.time()
  time    <- secs(t2, t1)
  if (!inherits(pout, "error")) {
    n_trial <- julia_eval(glue('p{routine}.trials'))
    success <- julia_eval(glue('p{routine}.convergence'))
  } else {
    error <- pout$message
  }

  # Visualise outputs
  if (FALSE) {
    # # Depth time series
    # arc <- args$.yobs$ModelObsDepthNormalTrunc
    # plot(arc$timestamp[1:600], arc$obs[1:600] * -1, type = "l")
    # points(pout$states$timestamp, pout$states$map_value * -1, pch = ".", col = scales::alpha("dimgrey", 0.5))
    # lines(arc$timestamp[1:600], arc$obs[1:600] * -1, lwd = 2)
    # # Depth 'error' (depth of individual below the seabed)
    # # > TO DO
    # # Animated maps 
    # m <- terra::rast(here_data("spatial", "bathy.tif"))
    # ani(.sim = sim, .map = m, .moorings = moorings, .start = 1L, .end = 600, .input = args, .output = pout)
    # # Interactive debugging
    # debug_acdc(.map = terra::rast(here_data("spatial", "bathy.tif")), .moorings = moorings, .step = 100, .input = args, .output = pout)
  }
  
  # Collect success statistics
  dout    <- data.table(id = sim$index, 
                        method = "particle", 
                        routine = routine, 
                        success = success, 
                        error = error, 
                        n_trial = n_trial,
                        time = time)
  
  # Write particles
  # * Particles currently exist in memory for smoothing
  # * Each particles dataset (for a one-month period) is ~314 MB 
  # * This is too big to write to file
  # if (success) {
  #  qs::qsave(pout, file.path(sim$folder_coord, pfile))
  # }
  
  # Write convergence
  qs::qsave(dout, file.path(sim$folder_coord, dfile))
  
  # Return success
  # * This controls whether or not we run smoothing
  success
}

#' Create a smoothing wrapper that runs the smoother for a given simulation
pf_smoother_wrapper <- function(sim) {
  
  # Define output files
  pfile   <- "coord-smo.qs"
  dfile   <- "data-smo.qs"
  
  # Run smoother
  t1      <- Sys.time()
  pout    <- pf_smoother_two_filter(.n_particle = 500L,
                                    .verbose = TRUE)
  t2 <- Sys.time()
  
  # Collect success statistics
  time <- secs(t2, t1)
  dout <- data.table(id = sim$index, 
                     method = "particle", 
                     routine = "smo", 
                     success = TRUE, 
                     error = NA_character_, 
                     time = time)
  
  # Write outputs
  # * Smoother outputs seem to be smaller than those from the filter (13 MB)
  # * This may be due to more efficient compression due to NAs (?)
  # * For now, this appears to be manageable. 
  # * (optional) TO DO: check file sizes across simulations
  qs::qsave(pout, file.path(sim$folder_coord, pfile))
  qs::qsave(dout, file.path(sim$folder_coord, dfile))
  nothing()
  
}