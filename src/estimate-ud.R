estimate_ud_spatstat <- function(sim, extract_coord = NULL, map, win, sigmas, plot = plot) {
  
  cat("\n\n\n---------------------------------------------------------------\n")
  msg("\n On row {sim$index}...", .envir = environment())
  
  # Define coordinates
  # sim <- iteration[1, ]
  stopifnot(rlang::has_name(sim, "file_coord"))
  coord <- qs::qread(sim$file_coord)
  if (!is.null(extract_coord)) {
    coord <- extract_coord(coord)
  }
  
  # Estimate maps using selected sigma methods
  shortcut <- list()
  for (sigma_index in seq_len(length(sigmas))) {
    
    #### Compute UD
    # UD list
    t1 <- Sys.time()
    ud_list <- map_dens(.map = map,
                        .owin = win, 
                        .coord = coord, 
                        .discretise = TRUE, 
                        .shortcut = shortcut,
                        sigma = sigmas[[sigma_index]],
                        .plot = plot,
                        .use_tryCatch = TRUE, 
                        .verbose = TRUE)
    # terra::plot(ud_list$ud)
    t2 <- Sys.time()
    # Convergence
    convergence <- !is.null(ud_list$ud) & !inherits(ud_list$ud, "error")
    # Time
    sim[, time := secs(t2, t1)]
    sim[, convergence := convergence]
    
    #### Write UD files to file
    # (optional) Write shortcut list
    folder_ud <- file.path(sim$folder, "ud", names(sigmas)[sigma_index])
    if (!is.null(ud_list$x) & !inherits(ud_list$x, "error") & 
        !is.null(ud_list$D) & !inherits(ud_list$D, "error")) {
      shortcut <- list(x = ud_list$x, D = ud_list$D)
      qs::qsave(shortcut, 
                file.path(folder_ud, "shortcut.qs"))
    }
    # Write convergence/timings
    qs::qsave(sim, file.path(folder_ud, "time.qs"))
    # Write UD
    if (convergence) {
      terra::writeRaster(ud_list$ud, file.path(folder_ud, "ud.tif"), overwrite = TRUE)
    } 
  }
  nothing()
}

estimate_ud_dbbmm <- function(sim, map, bbrast_ll) {

  
  cat("\n\n\n---------------------------------------------------------------\n")
  msg("\n On row {sim$index}...", .envir = environment())
  
  # Define coordinates
  # sim <- iteration[1, ]
  stopifnot(rlang::has_name(sim, "file_coord"))
  coord   <- qs::qread(sim$file_coord)
  success <- TRUE
  
  # Generate UD
  t1 <- Sys.time()
  dbb <- tryCatch(
    RSP::dynBBMM(input = coord,
                 base.raster = bbrast_ll,
                 UTM = "29"),
    error = function(e) e)
  t2 <- Sys.time()
  
  # Error handling
  if (inherits(dbb, "error")) {
    success <- FALSE
    message(dbb$message)
  }
  
  # Process raster
  if (success) {
    dbb <- dbb$dbbmm[[1]]
    dbb <- terra::rast(dbb)
    if (terra::nlyr(dbb) > 1L) {
      dbb <- spatNormalise(terra::app(dbb, "sum"))
    }
    # Resample RSP onto map grid for consistency
    ud <- terra::project(dbb, terra::crs(map))
    ud <- terra::resample(ud, map)
    ud <- terra::mask(ud, map)
    ud <- spatNormalise(ud)
  }

  # Save outputs
  if (success) {
    # TO DO
  }
  
}

lapply_estimate_ud_spatstat <- function(iteration, 
                                        extract_coord = NULL, 
                                        plot = FALSE, 
                                        cl = 2L) {
  
  tictoc::tic()
  on.exit(tictoc::toc(), add = TRUE)
  if (spatstat.geom::spatstat.options("npixel") < 500) {
    warning("Low pixel resolution!", immediate. = TRUE, call. = FALSE)
  }
  
  #### Load data for UD estimation
  # Read estimation window
  win <- qs::qread(dv::here_data("spatial", "win.qs"))
  # Define bandwidth estimators
  # * null     :    fast
  # * bw.ppl   : slow
  # * bw.diggle: slow
  sigmas <- list(null = NULL) 
                 # bw.ppl = spatstat.explore::bw.ppl,
                 # bw.diggle = spatstat.explore::bw.diggle)
  # Define UD function
  estimate_ud <- function(.xi, .chunkargs) {
    estimate_ud_spatstat(sim = .xi, 
                         extract_coord = extract_coord, 
                         map = .chunkargs, 
                         win = win, 
                         sigmas = sigmas, 
                         plot = plot)
  }
  
  # Build directories for UDs
  # {individual}/{mmyy}/{algorithm}/{parameter combination}/ud/{bandwidth}/
  utils.add::check_names(input = iteration, req = "folder")
  lapply(names(sigmas), function(sig) {
    dirs.create(file.path(iteration$folder, "ud", sig))
  })
  
  #### Estimate UDs for iterations in parallel 
  cl_lapply(.x = split(iteration, collapse::seq_row(iteration)), 
            .cl = cl,
            .chunk = TRUE, 
            .chunk_fun = function(...) {
              .chunkargs <- terra::rast(dv::here_data("spatial", "ud-grid.tif"))
              terra:::readAll(.chunkargs)
              .chunkargs
            },
            .fun = estimate_ud)
  nothing()
}

lapply_estimate_ud_dbbmm <- function() {
  
}