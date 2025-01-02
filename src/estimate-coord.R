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
               .detections = detections,
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
estimate_coord_patter <- function(sim, datasets, trial = FALSE) {
  
  #### Initialise
  coffee()
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
  # (For historical convenience, xinit is read below)
  # xinit       <- datasets$xinit_by_unit[[sim$unit_id]]
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
    warn(paste("> The timeline is", length(timeline), "steps long."))
  }
  
  #### Define initial locations
  # simulation: data/input/simulation/unit_id/xinit-fwd.qs, xinit-bwd.qs
  # real      : input/xinit/1.5/individual_id/month_id/xinit-fwd.qs, xinit-bwd.qs
  xinit_fwd <- qs::qread(file.path(sim$folder_xinit, "xinit-fwd.qs"))
  xinit_bwd <- qs::qread(file.path(sim$folder_xinit, "xinit-bwd.qs"))
  # Add behaviour
  # * For simulations, this is achieved in simulate-algorithms.R
  # * TO DO For real-world analyses, this could be achieved in prepare-xinits-*.R
  xinit_fwd[, behaviour := as.integer(behaviour[1])]
  xinit_bwd[, behaviour := as.integer(behaviour[length(behaviour) - 1L])]
  xinit_list <- list(forward  = xinit_fwd, 
                     backward = xinit_bwd)
  xinit_list <- list(forward  = xinit_fwd, 
                     backward = xinit_bwd)
  
  #### Define movement model
  state       <- state_flapper
  model_move  <- patter_ModelMove(sim)
  julia_assign("behaviour", behaviour)
  update_model_move_components()

  #### Define observation model(s)
  model_obs <- patter_ModelObs(sim        = sim, 
                               timeline   = timeline, 
                               detections = detections,
                               moorings   = moorings,
                               archival   = archival, 
                               xinit_list = xinit_list)
  # Fix ModelObsAcousticContainer elements
  # * All algorithms contain an 'ModelObsAcousticContainer' element
  # * Because for speed ModelObsCaptureContainer is also coded as ModelObsAcousticContainer
  yobs_fwd <- yobs_bwd <- model_obs
  yobs_fwd$ModelObsAcousticContainer <- model_obs$ModelObsAcousticContainer$forward
  yobs_bwd$ModelObsAcousticContainer <- model_obs$ModelObsAcousticContainer$backward
  stopifnot(!is.null(yobs_fwd$ModelObsAcousticContainer))
  stopifnot(!is.null(yobs_bwd$ModelObsAcousticContainer))
  
  #### List filter arguments
  # Set arguments to reduce computation time 
  args <- list(.timeline   = timeline,
               .state      = state,
               .xinit      = xinit_fwd,
               .yobs       = yobs_fwd,
               .model_move = model_move,
               .n_particle = sim$np, # patter_np(sim),
               .n_move     = 10000L,
               .n_record   = 3000L,  # 1500L for simulations, 2000L for real-world analysis
               .n_iter     = 1L,
               .direction  = "forward",
               .collect    = FALSE,
               .verbose    = TRUE)
  
  #### Implement forward filter
  cat("... (1) Implementing forward filter...\n")
  success <- pf_filter_wrapper(sim = sim, args = args)
  if (trial) {
    return(success)
  }
  
  #### Implement backward filter 
  if (success) {
    cat("\n... (2) Implementing backward filter...\n")
    
    #### Update behaviour
    
    # On forward run, behaviour at t - 1 is based on behaviour from t - 1 to t
    # * E.g., at t3 behaviour is defined based on behaviour from t3 to t4
    
    # On backward filter (and smoother), we move from t -> t-1 -> t-2 etc.
    # * We need to correct the behavioural vector accordingly
    # * I.e., at t3 behaviour is defined based on behaviour from t3 to t2
    
    # Time steps: 1,        2,        3,        4,        5
    # Forward   : 1 (1->2), 1 (2->3), 2 (3->4), 2 (4->5), 2 (5->)
    # Backward  : 2 (1->),  1 (2->1), 1 (3->2), 2 (4->3), 2 (5->4)
    # c(1, 1, 2, 2, 2)
    # dplyr::lag(c(1, 1, 2, 2, 2))
    
    behaviour    <- dplyr::lag(behaviour)
    behaviour[1] <- 2L
    julia_assign("behaviour", behaviour)
    update_model_move_components()
    
    #### Update args
    args$.xinit     <- xinit_bwd
    args$.yobs      <- yobs_bwd
    args$.direction <- "backward"
    
    #### Run filter
    success         <- pf_filter_wrapper(sim = sim, args = args)
  }

  #### Implement smoothing 
  if (success & sim$smooth) {
    cat("\n... (3) Implementing smoother...\n")
    success <- pf_smoother_wrapper(sim)
  }
  
  #### Clean up
  pfwd <- patter:::name_particles(.fun = "pf_filter", .direction = "forward")
  pbwd <- patter:::name_particles(.fun = "pf_filter", .direction = "backward")
  ptf  <- patter:::name_particles(.fun = "pf_smoother_two_filter")
  julia_command(glue('{pfwd} = nothing;'))
  julia_command(glue('{pbwd} = nothing;'))
  julia_command(glue('{ptf} = nothing;'))
  nothing()
}