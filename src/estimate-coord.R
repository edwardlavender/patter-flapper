#' Estimate COAs
estimate_coord_coa <- function(sim, map, datasets) {
  
  # Initialise
  cat_init(sim$index)
  
  # Define outputs
  pfile <- "coord.qs"
  dfile <- "data.qs"
  
  # Define data & parameters
  detections <- datasets$detections_by_unit[[sim$unit_id]]
  moorings   <- datasets$moorings
  delta_t    <- sim$delta_t
  
  # Run algorithm
  t1    <- Sys.time()
  pout  <- coa(.map = map,
               .acoustics = detections,
               .moorings = moorings,
               .delta_t = delta_t,
               .plot_weights = FALSE)
  t2    <- Sys.time()
  time  <- secs(t2, t1)
  
  # Collect success statistics
  dout <- data.table(id = sim$index, 
                     method = "coa", 
                     routine = "coa", 
                     success = TRUE,
                     error = NA_character_, 
                     ntrial = NA_integer_,
                     time = time)
  
  # Write outputs
  qs::qsave(pout, file.path(sim$folder_coord, pfile))
  qs::qsave(dout, file.path(sim$folder_coord, dfile))
  nothing()
  
}

#' Estimate RSPs
estimate_coord_rsp <- function(sim, map, datasets) {
  
  # Initialise
  cat_init(sim$index)
  
  # Define outputs
  pfile <- "coord.qs"
  dfile <- "data.qs"
  
  # Define data & parameters
  detections <- datasets$detections_by_unit[[sim$unit_id]]
  moorings   <- datasets$moorings
  act        <- as_actel(.map = map, 
                         .acoustics = detections, 
                         .moorings = moorings)
  tm         <- datasets$tm
  er.ad      <- sim$er.ad
  
  # Initialise algorithm
  error   <- NA_character_
  success <- TRUE
  args    <- list(input = act,
                  t.layer = tm,
                  coord.x = "Longitude", coord.y = "Latitude",
                  er.ad = er.ad)
  
  # Run algorithm
  t1 <- Sys.time()
  pout <- tryCatch(do.call(RSP::runRSP, args),
                   error = function(e) e)
  t2 <- Sys.time()
  
  # Handle errors
  if (inherits(pout, "error")) {
    success <- FALSE
    error   <- pout$message
    message(error)
  } else {
    time <- secs(t2, t1)
  }
  
  # Collect success statistics
  dout <- data.table(id = sim$index, 
                     method = "rsp", 
                     routine = "rsp", 
                     success = success, 
                     error = error, 
                     ntrial = NA_integer_,
                     time = time)
  
  # Write outputs
  if (success) {
    qs::qsave(pout, file.path(sim$folder_coord, pfile))
  }
  qs::qsave(dout, file.path(sim$folder_coord, dfile))
  nothing()
  
}

# Estimate particles
estimate_coord_patter <- function(sim, map, datasets) {
  
  #### Initialise
  cat_init(sim$index)
  
  #### Define outputs
  # See pf_filter_wrapper()
  # See pf_smoother_wrapper()
  
  #### Define inputs (data, parameters)
  # Parameters are defined by the patter_*() functions below
  detections  <- datasets$detections_by_unit[[sim$unit_id]]
  moorings    <- datasets$moorings
  archival    <- datasets$archival_by_unit[[sim$unit_id]]
  behaviour   <- datasets$behaviour_by_unit[[sim$unit_id]]
  # (optional) Trial implementation of depth data only when resting
  # archival <- archival[which(behaviour == 1L), ]
  # stopifnot(nrow(archival) > 0L)
  
  #### Define timeline
  timeline <- patter_timeline(sim$month_id)
  if (FALSE) {
    # Use a restricted timeline for testing
    warn("Using a restricted timeline for testing!")
    sel <- 1:10
    if (!is.null(detections)) {
      # Use the first N observations required to achieve at least two detections at different receivers
      # This is necessary so that $yobs$ModelObsAcousticContainer is not empty
      detections_timestamps <- lubridate::round_date(detections$timestamp, "2 mins")
      sel <- which(timeline %in% detections_timestamps[detections$receiver_id != detections$receiver_id[1]])
      sel <- 1:sel[1]
      timeline <- assemble_timeline(.datasets = plyr::compact(list(detections, archival)),
                                    .step = "2 mins",
                                    .trim = TRUE)
    }
    timeline <- timeline[sel]
  }
  warn(paste("> The timeline is", length(timeline), "steps long."))
  
  #### Define movement model
  state       <- "StateXY"
  model_move  <- patter_ModelMove(sim)
  julia_assign("behaviour", behaviour)
  JuliaCall::julia_command(simulate_step.ModelMoveFlapper)
  JuliaCall::julia_command(logpdf_step.ModelMoveFlapper)

  #### Define observation model(s)
  model_obs <- patter_ModelObs(sim = sim, 
                               timeline = timeline, 
                               detections = detections,
                               moorings = moorings,
                               archival = archival)
  
  yobs_fwd <- yobs_bwd <- model_obs
  if (rlang::has_name(model_obs, "ModelObsAcousticContainer")) {
    yobs_fwd$ModelObsAcousticContainer <- model_obs$ModelObsAcousticContainer$forward
    yobs_bwd$ModelObsAcousticContainer <- model_obs$ModelObsAcousticContainer$backward
  }
  
  #### List filter arguments
  # Set arguments to reduce computation time 
  args <- list(.timeline   = timeline,
               .state      = state,
               .xinit      = NULL,
               .yobs       = yobs_fwd,
               .model_move = model_move,
               .n_particle = sim$np, # patter_np(sim),
               .n_move     = 10000L,
               .n_record   = 1000L,
               .n_iter     = 1L,
               .direction  = "forward",
               .verbose    = TRUE)
  
  #### Implement particle filter
  # Forward filter 
  cat("... (1) Implementing forward filter...\n")
  success <- pf_filter_wrapper(sim = sim, args = args)
  # Backward filter 
  if (success) {
    cat("\n... (2) Implementing backward filter...\n")
    args$.yobs      <- yobs_bwd
    args$.direction <- "backward"
    success         <- pf_filter_wrapper(sim = sim, args = args)
  }

  #### Implement smoothing 
  if (success & sim$smooth) {
    cat("\n... (3) Implementing smoother...\n")
    success <- pf_smoother_wrapper(sim)
  }
  nothing()
  
}