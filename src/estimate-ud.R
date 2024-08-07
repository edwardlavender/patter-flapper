#' h bandwidth estimator

bw.h <- function(X) {
  # Code from adehabitatHR:::.kernelUDs with tweaks
  # & adapted for adehabitatHR
  sdxy <- sqrt(0.5 * (var(X$x) + var(X$y)))
  sdxy * (X$n^(-1/6))
}

if (FALSE) {
  library(testthat)
  library(adehabitatHR)
  data(puechabonsp)
  loc      <- puechabonsp$relocs
  ud       <- kernelUD(loc[, 1])
  h_expect <- ud$Brock@h$h
  coord    <- loc[loc$Name == "Brock", ]@coords
  X        <- data.frame(x = coord[, 1], y = coord[, 2])
  sdxy     <- sqrt(0.5 * (var(X$x) + var(X$y)))
  h_calc   <- sdxy * (nrow(X)^(-1/6))
  expect_equal(h_calc, h_expect)
  # > TRUE
}

#' Estimate UDs using spatstat
estimate_ud_spatstat <- function(sim, extract_coord = NULL, map, win, sigmas, plot) {
  
  # Initialise
  cat_init(sim$index)
  
  # Define outputs
  pfile <- "ud.tif"
  dfile <- "data.qs"
  
  # Define data 
  stopifnot(rlang::has_name(sim, "file_coord"))
  if (!file.exists(sim$file_coord)) {
    message("Coordinate file does not exist.")
    return(nothing())
  }
  coord <- qs::qread(sim$file_coord)
  if (!is.null(extract_coord)) {
    coord <- extract_coord(coord)
  }
  
  # Run algorithm 
  # (Estimate UDs for each sigma parameter)
  shortcut <- list()
  for (sigma_index in seq_len(length(sigmas))) {
    
    # Define parameters
    sigma       <- sigmas[[sigma_index]]
    sigma_label <- names(sigmas)[sigma_index]
    
    # Initialise algorithm
    error   <- NA_character_
    success <- TRUE
    
    # Estimate UD
    t1 <- Sys.time()
    ud <- tryCatch(
      {
        map_dens(.map = map,
                 .owin = win, 
                 .coord = coord, 
                 .discretise = TRUE, 
                 .shortcut = shortcut,
                 sigma = sigma,
                 .plot = plot,
                 .use_tryCatch = FALSE, 
                 .verbose = TRUE) 
      }, 
      error = function(e) e)
    t2     <- Sys.time()
    
    # Error handling
    if (inherits(ud, "error")) {
      success <- FALSE
      error   <- ud$message
      message(error)
    } else {
      time <- secs(t2, t1)
    }
    
    # Collect success statistics
    dout <- data.table(id = sim$index, 
                       method = "spatstat", 
                       routine = sigma, 
                       success = success,
                       error = error, 
                       convergence = success, 
                       ntrial = NA_integer_,
                       time = time)
    
    # Write outputs
    folder_ud <- file.path(sim$folder_ud, "spatstat", sigma_label)
    if (success) {
      shortcut <- list(x = ud$x, D = ud$D)
      qs::qsave(shortcut, file.path(folder_ud, "shortcut.qs"))
      terra::writeRaster(ud$ud, 
                         file.path(folder_ud, pfile), 
                         overwrite = TRUE)
    }
    qs::qsave(sim, file.path(folder_ud, dfile))
    
  }
  
  nothing()
  
}

#' Estimate DBBMMs
estimate_ud_dbbmm <- function(sim, map, bbrast_ll, plot) {
  
  # Initialise
  cat_init(sim$index)
  
  # Define outputs
  pfile <- "ud.tif"
  dfile <- "data.qs"
  
  # Define data 
  stopifnot(rlang::has_name(sim, "file_coord"))
  if (!file.exists(sim$file_coord)) {
    message("Coordinate file does not exist.")
    return(nothing())
  }
  coord   <- qs::qread(sim$file_coord)

  # Initialise algorithm
  error   <- NA_character_
  success <- TRUE
  
  # Run algorithm
  t1 <- Sys.time()
  ud <- tryCatch(
    {
    # Estimate UD
    dbb <- RSP::dynBBMM(input = coord,
                        base.raster = bbrast_ll,
                        UTM = "29")
    # Process UD
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
    if (plot) {
      terra::plot(ud)
    }
    ud
  }, error = function(e) e)
  t2 <- Sys.time()
  
  # Error handling
  if (inherits(ud, "error")) {
    success <- FALSE
    error   <- ud$message
    message(error)
  } else {
    time <- secs(t2, t1)
  }
  
  # Collect success statistics
  dout <- data.table(id = sim$index, 
                     method = "rsp", 
                     routine = "dbbmm", 
                     success = success,
                     error = error, 
                     convergence = success, 
                     ntrial = NA_integer_,
                     time = time)
  
  # Write outputs
  folder_ud <- file.path(sim$folder_ud, "dbbmm")
  if (success) {
    terra::writeRaster(ud, 
                       file.path(folder_ud, pfile), 
                       overwrite = TRUE)
  }
  qs::qsave(sim, file.path(folder_ud, dfile))
  nothing()
  
}