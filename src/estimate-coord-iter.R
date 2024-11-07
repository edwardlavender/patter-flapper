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
lapply_estimate_coord_patter <- function(iteration, datasets, trial = FALSE) {
  
  tictoc::tic()
  on.exit(tictoc::toc(), add = TRUE)
  
  # Estimate coordinates via estimate_coord_patter()
  iteration_ls <- split(iteration, collapse::seq_row(iteration))
  cl_lapply(iteration_ls, function(sim) {
    estimate_coord_patter(sim = sim, 
                          datasets = datasets, 
                          trial = trial)
  })
  
  nothing()
  
}
