estimate_uds <- function(iteration) {
  
  #### Load data for UD estimation
  # Read estimation window
  win <- qs::qread(here_data("spatial", "win.qs"))
  # Define bandwidth estimators
  sigmas <- list(bw.ppl = spatstat.explore::bw.ppl,
                 bw.diggle = spatstat.explore::bw.diggle)
  
  # Build directories for UDs
  # {individual}/{mmyy}/{algorithm}/{parameter combination}/ud/{bandwidth}/
  lapply(names(sigmas), function(sig) {
    dirs.create(file.path(iteration$folder, "ud", sig))
  })
  
  #### Estimate UDs for iterations in parallel 
  cl_lapply(.x = iteration_ls[1:2], 
            .cl = 2L,
            .chunk = TRUE, 
            .chunk_fun = function(...) {
              .chunkargs <- terra::rast(here_data("spatial", "ud-grid.tif"))
              terra:::readAll(.chunkargs)
              .chunkargs
            },
            .fun = function(sim, .chunkargs) {
              
              # Define coordinates
              # sim <- iteration[1, ]
              coord <- qs::qread(file.path(sim$folder, "coord.qs"))
              
              # Estimate maps using selected sigma methods
              shortcut <- list()
              for (sigma_index in seq_len(length(sigmas))) {
                
                #### Compute UD
                # UD list
                t1 <- Sys.time()
                ud_list <- map_dens(.map = .chunkargs,
                                    .owin = win, 
                                    .coord = coord, 
                                    .discretise = TRUE, 
                                    .shortcut = shortcut,
                                    sigma = sigmas[[sigma_index]],
                                    .use_tryCatch = FALSE, 
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
                
                nothing()
                
              }
              
              nothing()
              
            })
  
}