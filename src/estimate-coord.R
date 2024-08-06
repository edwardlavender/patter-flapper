estimate_coord_coa <- function(sim, map, datasets) {
  
  # Estimate COAs
  t1 <- Sys.time()
  coord <- coa(.map = map,
               .acoustics = datasets$detections_by_unit[[sim$unit_id]],
               .moorings = datasets$moorings,
               .delta_t = sim$delta_t,
               .plot_weights = FALSE)
  t2 <- Sys.time()
  
  # Record convergence & timing
  # print(head(coord))
  time <- sim[, time := secs(t2, t1)]
  
  # Save outputs
  qs::qsave(time, file.path(sim$folder, "data.qs"))
  qs::qsave(coord, file.path(sim$folder, "coord.qs"))
  nothing()
  
}

estimate_coord_rsp <- function() {
  
  
  
}

estimate_coord_patter <- function(sim, map, datasets) {
  
  #### Read data
  # sim <- iteration[81, ]
  cat("\n\n\n---------------------------------------------------------------\n")
  msg("\n On row {sim$index}...", .envir = environment())
  detections  <- datasets$detections_by_unit[[sim$unit_id]]
  archival    <- datasets$archival_by_unit[[sim$unit_id]]
  moorings    <- datasets$moorings

  #### Particle filter arguments 
  # Timeline
  timeline <- patter_timeline(sim$month_id)
  if (TRUE) {
    # Use a restricted timeline for testing
    warn("Using a restricted timeline for testing!")
    timeline <- assemble_timeline(.datasets = plyr::compact(list(detections, archival)),
                                  .step = "2 mins",
                                  .trim = TRUE)[1:10L]
  }
  # range(timeline)
  # Movement model
  state       <- "StateXY"
  model_move  <- patter_ModelMove(sim)
  # Observation model(s)
  model_obs <- patter_ModelObs_forward(sim = sim, 
                                       timeline = timeline, 
                                       detections = detections,
                                       moorings = moorings,
                                       archival = archival)
  # List arguments
  args <- list(.map = map,
               .timeline = timeline,
               .state = state,
               .xinit = NULL,
               .xinit_pars = list(mobility = sim$mobility),
               .yobs = model_obs$yobs,
               .model_obs = model_obs$model_obs,
               .model_move = model_move,
               .n_particle = patter_np(sim),
               .direction = "forward",
               .verbose = TRUE)
  
  #### Implement particle filter
  # Forward filter 
  cat("... (1) Implementing forward filter...\n")
  success <- pf_filter_wrapper(sim = sim, args = args)
  # Backward filter 
  if (success) {
    cat("\n ... (2) Implementing backward filter...\n")
    args$.yobs      <- patter_ModelObs_backward(sim, model_obs)$yobs
    args$.direction <- "backward"
    success         <- pf_filter_wrapper(sim = sim, args = args)
  }

  #### Implement smoothing 
  if (success) {
    cat("\n... (3) Implementing smoother...\n")
    success <- pf_smoother_wrapper(sim)
  }

  #### Return outputs
  success
  
}

lapply_estimate_coord_coa <- function(iteration, datasets) {
  
  tictoc::tic()
  on.exit(tictoc::toc(), add = TRUE)
  
  # Read grid
  map <- terra::rast(dv::here_data("spatial", "ud-grid.tif"))
  terra:::readAll(map)
  
  # Estimate COAs iteratively
  iteration_ls <- split(iteration, collapse::seq_row(iteration))
  cl_lapply(iteration_ls, function(sim) {
    estimate_coord_coa(sim = sim, 
                       map = map,
                       datasets = datasets)
  })
  
  nothing()
  
}


lapply_estimate_coord_patter <- function(iteration, datasets) {
  
  tictoc::tic()
  on.exit(tictoc::toc(), add = TRUE)
  
  # Read grid
  map <- terra::rast(dv::here_data("spatial", "bathy.tif"))
  terra:::readAll(map)
  
  # Estimate coordinates iteratively
  iteration_ls <- split(iteration, collapse::seq_row(iteration))
  cl_lapply(iteration_ls, function(sim) {
    estimate_coord_patter(sim = sim, 
                          map = map,
                          datasets = datasets)
  })
  
  nothing()
  
}

