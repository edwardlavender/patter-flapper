if (!patter:::os_linux()) {
 
  #' Estimate COAs iteratively
  lapply_estimate_coord_coa <- function(iteration, datasets) {
    
    tictoc::tic()
    on.exit(tictoc::toc(), add = TRUE)
    
    # Read grid
    map <- terra::rast(dv::here_data("spatial", "ud-grid.tif"))
    terra:::readAll(map)
    
    # Estimate COAs via estimate_coord_coa()
    iteration_ls <- split(iteration, collapse::seq_row(iteration))
    cl_lapply(iteration_ls, function(sim) {
      estimate_coord_coa(sim = sim, 
                         map = map,
                         datasets = datasets)
    })
    
    nothing()
    
  }
   
  #' Estimate RSPs iteratively
  lapply_estimate_coord_rsp <- function(iteration, datasets) {
    
    tictoc::tic()
    on.exit(tictoc::toc(), add = TRUE)
    
    # Read grid 
    map <- terra::rast(dv::here_data("spatial", "bathy.tif"))
    terra:::readAll(map)
    
    # Read tm
    datasets$tm <- qs::qread(dv::here_data("spatial", "tm.qs"))
    
    # Estimate coordinates via estimate_coord_rsp()
    iteration_ls <- split(iteration, collapse::seq_row(iteration))
    cl_lapply(iteration_ls, function(sim) {
      estimate_coord_rsp(sim = sim, 
                         map = map,
                         datasets = datasets)
    })
    
    nothing()
    
  }
  
}

#' Estimate particles iteratively
lapply_estimate_coord_patter <- function(iteration, 
                                         datasets, 
                                         trial = FALSE, 
                                         log.folder = NULL, 
                                         log.txt = NULL) {
  
  # (optional) Open sink
  log.txt <- sink_open(log.folder = log.folder, log.txt = log.txt)
  on.exit(sink_close(log.txt), add = TRUE)
  
  # Print time 
  print(Sys.time())
  tictoc::tic()
  on.exit(tictoc::toc(), add = TRUE)
  
  # Clean up coffee-break directory 
  unlink(list.files(here_data("coffee"), full.names = TRUE))
  
  # Estimate coordinates via estimate_coord_patter()
  iteration_ls <- split(iteration, collapse::seq_row(iteration))
  success <- cl_lapply(iteration_ls, function(sim) {
    estimate_coord_patter(sim = sim, 
                          datasets = datasets, 
                          trial = trial)
  })
  
  # Close function
  print(Sys.time())
  if (trial) {
    return(unlist(success))
  }
  nothing()
  
}
