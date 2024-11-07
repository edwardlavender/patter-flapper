if (!patter:::os_linux()) {
  
#' Estimate UDs via spatstat iteratively 
lapply_estimate_ud_spatstat <- function(iteration, 
                                        extract_coord = NULL, 
                                        plot = FALSE, 
                                        cl = 2L) {
  
  #### Initialise
  tictoc::tic()
  on.exit(tictoc::toc(), add = TRUE)
  if (spatstat.geom::spatstat.options("npixel") < 500) {
    warning("Low pixel resolution!", immediate. = TRUE, call. = FALSE)
  }
  
  #### Define data & parameters
  # Read estimation window
  win <- qs::qread(dv::here_data("spatial", "win.qs"))
  # Define bandwidth estimators
  # * h, null (fast), bw.ppl (slow), bw.diggle (slow)
  sigmas <- list(h = bw.h) 
  # spatstat/{sigma} folders should exist 
  # * These are defined in R/prepare-runs.R
  sapply(names(sigmas), function(sigma_label) {
    stopifnot(dir.exists(file.path(iteration$folder_ud, "spatstat", sigma_label)))
  })
  # Define UD function
  estimate_ud <- function(.xi, .chunkargs) {
    estimate_ud_spatstat(sim = .xi, 
                         extract_coord = extract_coord, 
                         map = .chunkargs, 
                         win = win, 
                         sigmas = sigmas, 
                         plot = plot)
  }
  
  #### Estimate UDs iteratively
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

#' Estimate DBBMMs iteratively 
lapply_estimate_ud_dbbmm <- function(iteration, 
                                     plot = FALSE, 
                                     cl = 2L) {
  
  #### Initialise
  tictoc::tic()
  on.exit(tictoc::toc(), add = TRUE)
  
  #### Define data & parameters
  estimate_ud <- function(.xi, .chunkargs) {
    estimate_ud_dbbmm(sim = .xi, 
                      map = .chunkargs$map, 
                      bbrast_ll = .chunkargs$bbrast_ll,
                      plot = plot)
  }
  
  #### Estimate UDs iteratively
  cl_lapply(.x = split(iteration, collapse::seq_row(iteration)), 
            .cl = cl,
            .chunk = TRUE, 
            .chunk_fun = function(...) {
              map       <- terra::rast(dv::here_data("spatial", "ud-grid.tif"))
              bbrast_ll <- terra::rast(dv::here_data("spatial", "bbrast_ll.tif"))
              terra:::readAll(map)
              terra:::readAll(bbrast_ll)
              list(map = map, bbrast_ll = bbrast_ll)
            },
            .fun = estimate_ud)
  
  nothing()
  
}

}
