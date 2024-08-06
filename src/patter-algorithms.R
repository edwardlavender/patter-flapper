patter_timeline <- function(mmyy) {
  period <- mmyyrng(mmyy)
  seq(period[1],  period[2], by = "2 mins")
}

patter_ModelMove <- function(sim) {
  move_xy(dbn_length =
            glue::glue("truncated(Gamma({sim$shape},
                                    {sim$scale}),
                              upper = {sim$mobility})"),
          dbn_angle = "Uniform(-pi, pi)")
}

patter_ModelObs_forward <- function(sim, timeline, detections, moorings, archival) {
  
  #### Acoustics datasets
  if (sim$dataset %in% c("ac", "acdc")) {
    # Detection model parameters
    receiver_alpha <- receiver_beta <- receiver_gamma <- NULL
    moorings[, receiver_alpha := sim$receiver_alpha]
    moorings[, receiver_beta := sim$receiver_beta]
    moorings[, receiver_gamma := sim$receiver_gamma]
    # Assemble acoustic data
    acoustics <- assemble_acoustics(.timeline = timeline,
                                    .acoustics = detections, 
                                    .moorings = moorings)
    # Assemble container data
    containers <- assemble_acoustics_containers(.acoustics = acoustics,
                                                .direction = "forward",
                                                .mobility = sim$mobility)
  }
  
  #### Archival datasets
  if (sim$dataset %in% c("dc", "acdc")) {
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
  }
  
  #### Return list of algorithm-specific ModelObs strings & corresponding datasets
  if (sim$dataset == "ac") {
    out <- list(model_obs = c("ModelObsAcousticLogisTrunc", 
                              "ModelObsAcousticContainer"), 
                yobs = list(ModelObsAcousticLogisTrunc = acoustics, 
                            ModelObsAcousticContainer = containers))
  } else if (sim$dataset == "dc") {
    out <- list(model_obs = "ModelObsDepthNormalTrunc", 
                yobs = list(ModelObsDepthNormalTrunc = archival))
  } else if (sim$dataset == "acdc") {
    out <- list(model_obs = c("ModelObsAcousticLogisTrunc", 
                              "ModelObsAcousticContainer", 
                              "ModelObsDepthNormalTrunc"), 
                yobs = list(ModelObsAcousticLogisTrunc = acoustics,
                            ModelObsAcousticContainer = containers, 
                            ModelObsDepthNormalTrunc = archival))
  }
  
  out
  
}

patter_ModelObs_backward <- function(sim, model_obs) {
  if (sim$dataset %in% c("ac", "acdc")) {
    # Update container dataset (which is direction-dependent)
    model_obs$yobs$ModelObsAcousticContainer <- 
      assemble_acoustics_containers(.acoustics = model_obs$yobs$ModelObsAcousticLogisTrunc,
                                    .direction = "backward",
                                    .mobility = sim$mobility)
  }
  model_obs
} 

patter_np <- function(sim) {
  if (sim$dataset == "ac") {
    np <- 1e4L
  } else if (sim$dataset == "dc") {
    np <- 1e5L
  } else if (sim$dataset == "acdc") {
    np <- 1e5L
  }
  np
}

pf_filter_init <- function(.map,
                           .timeline,
                           .state = "StateXY", .xinit = NULL, .xinit_pars = list(),
                           .yobs = list(),
                           .model_obs,
                           .model_move = move_xy(),
                           .n_move = 1e5L,
                           .n_particle = 1000L,
                           .n_resample = as.numeric(.n_particle),
                           .n_record = 1e3L,
                           .direction = c("forward", "backward"),
                           .verbose = getOption("patter.verbose")) {
  
  #### Initiate
  cats <- cat_setup(.fun = "pf_filter_init", .verbose = .verbose)
  on.exit(eval(cats$exit, envir = cats$envir), add = TRUE)
  
  #### Check user inputs
  cats$cat(paste0("... ", call_time(Sys.time(), "%H:%M:%S"), ": Checking user inputs..."))
  check_inherits(.state, "character")
  check_inherits(.model_move, "character")
  .direction <- match.arg(.direction)
  tzs <- c(tz(.timeline), sapply(.yobs, \(d) tz(d$timestamp)))
  if (length(unique(tzs)) != 1L) {
    abort("There is a mismatch between the time zones of `.timeline` and/or `.yobs` `timestamp`s ({str_items(tzs, .quo = '\"')}).",
          .envir = environment())
  }
  
  #### Set initial state
  cats$cat(paste0("... ", call_time(Sys.time(), "%H:%M:%S"), ": Setting initial states..."))
  set_timeline(.timeline)
  .xinit <- sim_states_init(.map = .map,
                            .timeline = .timeline,
                            .direction = .direction,
                            .datasets = .yobs,
                            .models = .model_obs,
                            .pars = .xinit_pars,
                            .state = .state,
                            .xinit = .xinit,
                            .n = .n_particle)
  set_states_init(.xinit = .xinit, .state = .state)
  
  #### Set filter arguments
  cats$cat(paste0("... ", call_time(Sys.time(), "%H:%M:%S"), ": Setting observations..."))
  set_yobs_via_datasets(.datasets = .yobs, .model_obs = .model_obs)
  set_model_move(.model_move)
  nothing()
  
}

pf_filter_run <- function(args) {
  # Run filter
  pf_obj <- set_pf_filter(.n_move = 100000L,
                          .n_resample = as.numeric(args$.n_particle),
                          .n_record = 1000L,
                          .direction = args$.direction)
  # Get particles in R
  pf_particles(.xinit = NULL, .pf_obj = pf_obj)
}

pf_filter_wrapper <- function(sim, args) {
  
  # Define output files (fwd, bwd, smo) 
  if (args$.direction == "forward") {
    routine <- "fwd"
    pfile   <- "coord-fwd.qs"
    dfile   <- "data-fwd.qs"
  } else if (args$.direction == "backward") {
    routine <- "bwd"
    pfile   <- "coord-bwd.qs"
    dfile   <- "data-bwd.qs"
  } else {
    stop("args$.direction is missing!")
  }
  
  # Initialise variables
  success     <- FALSE
  error       <- NA_character_
  convergence <- FALSE
  time        <- NA_real_
  
  # Initialise filter
  # * We initialise the filter once, for speed
  cat("... ... (a) Initialising filter...\n\n")
  init <- tryCatch(do.call(pf_filter_init, args, quote = TRUE),
           error = function(e) e)
  if (inherits("init", "error")) {
    error <- init$message
    message(error)
  }
  
  # Iteratively run filter until error/convergence/n_trials is reached
  cat("\n... ... (b) Iteratively implementing filter...\n")
  if (!inherits(init, "error")) {
    
    n_trial     <- 1L
    n_trial_req <- NA_integer_
    while (n_trial < 10L) {
      
      cat(paste0("... ... ... On trial ", n_trial, "...\n"))
      
      # Implement the filter
      t1   <- Sys.time()
      pout <- tryCatch(pf_filter_run(args), error = function(e) e)
      t2   <- Sys.time()
      # Handle errors
      if (inherits(pout, "error")) {
        error   <- pout$message
        message(error)
        n_trial <- Inf
      } else {
        # If convergence failure, try again
        if (!pout$convergence) {
          n_trial <- n_trial + 1L
          # Otherwise, record time & close 
        } else {
          time        <- secs(t2, t1)
          n_trial_req <- n_trial
          n_trial     <- Inf
        }
      }
    }
    
  }
  
  # Collect success statistics
  success <- !inherits(pout, "error") && pout$convergence
  dout    <- data.table(routine = routine, 
                        success = success, 
                        error = error, 
                        convergence = convergence, 
                        n_trial = n_trial_req,
                        time = time)
  
  # Write outputs to file
  if (success) {
    qs::qsave(pout, file.path(sim$folder, pfile))
  }
  qs::qsave(dout, file.path(sim$folder, dfile))
  
  # Return success
  success
}

pf_smoother_wrapper <- function(sim) {
  
  # Define output files
  pfile   <- "coord-smo.qs"
  dfile   <- "data-smo.qs"
  
  # Run smoother
  t1      <- Sys.time()
  pout <- pf_smoother_two_filter(.n_particle = 1000L,
                                 .verbose = FALSE)
  t2 <- Sys.time()
  
  # Collect success statistics
  time <- secs(t2, t1)
  dout <- data.table(routine = "smo", 
                     success = TRUE, 
                     error = NA_character_, 
                     convergence = TRUE, 
                     time = time)
  
  # Write outputs to file
  qs::qsave(pout, file.path(sim$folder, pfile))
  qs::qsave(dout, file.path(sim$folder, dfile))
  TRUE
  
}