###########################
###########################
#### run-patter-reanalysis.R

#### Aims
# 1) Try running ACDC algorithm with high resolution data
# ... to verify patterns in initial analysis

#### Prerequisites
# 1) Previous scripts


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
dv::src()

#### Load data 
if (!patter:::os_linux()) {
  map <- terra::rast(here_data("spatial", "bathy-5m.tif"))
}
pars              <- qs::qread(here_data("input", "pars.qs"))
iteration         <- qs::qread(here_data("input", "iteration", "patter.qs"))
moorings          <- qs::qread(here_data("input", "mefs", "moorings.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))
archival_by_unit  <- qs::qread(here_data("input", "archival_by_unit.qs"))
behaviour_by_unit <- qs::qread(here_data("input", "behaviour_by_unit.qs"))


###########################
###########################
#### Run algorithm 

#### Julia Set up
# Connect to Julia
if (patter:::os_linux()) {
  stopifnot(!any(c("terra", "sf") %in% sort(loadedNamespaces())))
}
julia_connect()
set_seed()
# Set map 
set_map(here_data("spatial", "bathy-5m.tif"))
# Set vmap & movement model components
batch <- pars$pmovement$mobility[1]
if (FALSE) {
  # Compute vmap (on MacOS): 1.7 mins
  tic()
  vmap <- spatVmap(.map = map, .mobility = batch, .plot = TRUE)
  terra::writeRaster(vmap, 
                     here_data("spatial", glue("vmap-5m-{batch}.tif")), 
                     overwrite = TRUE)
  toc()
} 
set_vmap(.vmap = here_data("spatial", glue("vmap-5m-{batch}.tif")))
set_model_move_components()

#### Prepare iterations/datasets
# Define iteration
iteration <- 
  iteration |>
  filter(dataset == "acdc" & sensitivity == "best") |> 
  mutate(
    # Update depth observation model for high-res bathymetry
    depth_sigma = 10, 
    depth_deep_eps = 40,
    # Update particles
    np = 5e6L,
    # Update folders
    folder_xinit = file.path("data", "input", "xinit", "1.5", individual_id, month_id),
    folder_coord = file.path("data", "output", "reanalysis", individual_id, month_id, "patter", "acdc", "1", "coord"), 
    file_coord   = file.path(folder_coord, "coord-smo.qs"),
    folder_ud    = file.path("data", "output", "reanalysis", individual_id, month_id, "patter", "acdc", "1", "ud")) |> 
  as.data.table()
# Define datasets
datasets <- list(detections_by_unit = acoustics_by_unit, 
                 moorings           = moorings,
                 archival_by_unit   = archival_by_unit, 
                 behaviour_by_unit  = behaviour_by_unit)
# Build dirs
if (FALSE) {
  unlink(iteration$folder_coord, recursive = TRUE)
  unlink(iteration$folder_ud, recursive = TRUE)
}
dirs.create(iteration$folder_coord)
# dirs.create(iteration$folder_ud)
# dirs.create(file.path(iteration$folder_ud, "spatstat", "h"))
dir.create(here_data("output", "log", "real", "reanalysis"), recursive = TRUE)

#### Estimate coordinates
if (TRUE) {
  
  #### Select rows
  nrow(iteration)
  rows <- 1:6     # 1
  # rows <- 7:12    # 2
  # rows <- 13:18   # 3
  # rows <- 19:24   # 4
  # rows <- 25:30   # 5
  # rows <- 31:36   # 6
  # rows <- 37:42   # 7 
  # rows <- 43:48   # 8
  iteration <- iteration[rows, ]
  
  #### Estimate coordinates
  gc()
  nrow(iteration)
  lapply_estimate_coord_patter(iteration  = iteration,
                               datasets   = datasets, 
                               trial      = FALSE, 
                               log.folder = here_data("output", "log", "real", "reanalysis"), 
                               log.txt    = glue("log-{min(rows)}.txt"))
  
}

#### Convergence 
# TO DO


#### End of code. 
###########################
###########################
