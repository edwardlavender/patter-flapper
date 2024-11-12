if (!patter:::os_linux()) {

#' Estimate residency from coordinates (particles)
lapply_estimate_residency_coord <- 
  function(files, extract_coord = NULL, cl = NULL) {
    
    # Set up
    tic()
    on.exit(toc(), add = TRUE)
    
    # Iteratively compute residency for files
    cl_lapply(
      .x = files, .cl = cl,.chunk = TRUE,
      .chunk_fun = function(file) {
        map <- terra::rast(here_data("spatial", "ud-grid.tif"))
        mpa <- qreadvect(here_data("spatial", "mpa.qs"))
        list(map = map, mpa = mpa)
      },
      .fun = function(file, .chunkargs) {
        
        # Set up
        # file = here_data("input", "simulation", "1", "coord.qs")
        map <- .chunkargs$map
        mpa <- .chunkargs$mpa
        
        # Handle convergence failures
        if (!file.exists(file)) {
          return(data.table(file = file, 
                            zone = c("closed", "open", "total"), 
                            time = NA_real_))
        }
        
        # Read coordinates
        coord <- qs::qread(file)
        if (!is.null(extract_coord)) {
          coord <- extract_coord(coord)
        }
        
        # Compute residency in open/closed regions (~10 s for patter)
        # tic()
        residency <- 
          coord |>
          # Discretise coordinates (for speed)
          select("x", "y") |>
          mutate(cell_id = terra::cellFromXY(map, cbind(x, y))) |> 
          group_by(cell_id) |>
          summarise(n = n()) |> 
          ungroup() |>
          # Determine if coordinates fall in open/closed areas
          mutate(x = terra::xFromCell(map, cell_id), 
                 y = terra::yFromCell(map, cell_id), 
                 open = terra::extract(mpa, cbind(x, y))$open, 
                 open = as.character(open)) |>
          # Sum up the number of coordinates in each area
          group_by(open) |> 
          summarise(time = sum(n)) |>
          ungroup() |> 
          # Compute proportion of time in each area
          mutate(time = time / sum(time)) |>
          filter(!is.na(open)) |>
          as.data.table()
        # toc()
        
        # Residency in MPA as a whole
        # sum(residency$time)
        
        # Return data.table with residency metrics
        data.table(file = file, 
                   zone = c(residency$open, "total"), 
                   time = c(residency$time, sum(residency$time)))
        
      }) |> rbindlist()
    
  }

#' Estimate residency from UDs (e.g., COA, RSPs)
lapply_estimate_residency_ud <- function(files) {
  
  # Setup 
  tic()
  on.exit(toc(), add = TRUE)
  stopifnot(all(file.exists(files)))
  mpa    <- qreadvect(here_data("spatial", "mpa.qs"))
  
  # Iteratively compute residency for files
  cl_lapply(files, function(file) {
    
    # Read UD
    # file <- here_data("input", "simulation", "1", "ud.tif")
    ud   <- terra::rast(file)
    names(ud) <- "lyr.1"
    
    # Compute residency in open/closed regions
    residency <- 
      terra::zonal(ud, mpa, "sum", na.rm = TRUE) |> 
      as.data.table() |> 
      mutate(id = mpa$id, 
             open = as.character(mpa$open)) |>
      select(id, open, time = "lyr.1") |> 
      group_by(open) |> 
      # Sum up probability densities in open/closed regions
      summarise(time = sum(time)) |>
      ungroup() |> 
      as.data.table()
    
    # Residency in MPA as a whole
    # sum(residency$time)
    
    # Return data.table with residency metrics
    data.table(file = file, 
               zone = c(residency$open, "total"), 
               time = c(residency$time, sum(residency$time)))
    
  }) |> rbindlist()
  
}

}