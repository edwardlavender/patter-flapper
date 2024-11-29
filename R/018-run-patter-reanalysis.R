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
dirs.create(iteration$folder_ud)
dirs.create(file.path(iteration$folder_ud, "spatstat", "h"))
dir.create(here_data("output", "log", "real", "reanalysis"), recursive = TRUE)

#### Estimate coordinates
if (FALSE) {
  
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

#### Patter error check 
# (Code copied from simulate-algorithms.R)
# Check for errors on forward filter 
sapply(split(iteration, seq_row(iteration)), function(d) {
  qs::qread(file.path(d$folder_coord, "data-fwd.qs"))$error
}) |> unlist() |> unique()
# Check for errors on backward filter 
sapply(split(iteration, seq_row(iteration)), function(d) {
  file <- file.path(d$folder_coord, "data-bwd.qs")
  if (file.exists(file)) {
    qs::qread(file)$error
  }
}) |> unlist() |> unique()
# Check for errors on smoother
sapply(split(iteration, seq_row(iteration)), function(d) {
  file <- file.path(d$folder_coord, "data-smo.qs")
  if (file.exists(file)) {
    qs::qread(file)$error
  }
}) |> unlist() |> unique()

#### Patter convergence 
iteration[, convergence := file.exists(file.path(folder_coord, "coord-smo.qs"))]
table(iteration$convergence)

#### Estimate UDs
if (FALSE && !patter:::os_linux()) {
  
  #### Estimate UDs
  # Time trial 
  lapply_estimate_ud_spatstat(iteration     = iteration[1, ], 
                              extract_coord = function(s) s$states,
                              cl            = NULL, 
                              plot          = FALSE)
  # Implementation (~1 min)
  lapply_estimate_ud_spatstat(iteration     = iteration, 
                              extract_coord = function(s) s$states,
                              cl            = 10L, 
                              plot          = FALSE)
  
}

#### Mapping
iteration |> 
  mutate(mapfile = file.path(folder_ud, "spatstat", "h", "ud.tif")) |>
  select(row = individual_id, column = month_id, mapfile) |> 
  ggplot_maps(png_args = list(filename = here_fig("analysis", "map-reanalysis.png"), 
                              height = 10, width = 9, units = "in", res = 800))

#### Residency
iteration_res <- copy(iteration)
iteration_res[, file := file.path(iteration$folder_ud, "spatstat", "h", "ud.tif")]
iteration_res[, file_exists := file.exists(file)]
residency <- lapply_estimate_residency_ud(files = iteration_res$file[iteration_res$file_exists])
# Write output
residency <- 
  left_join(iteration_res, residency, by = "file") |> 
  mutate(algorithm = dataset) |>
  select(individual_id, month_id, unit_id, algorithm, sensitivity, zone, time) |> 
  arrange(individual_id, month_id, unit_id, algorithm, sensitivity, zone) |>
  as.data.table()
qs::qsave(residency, here_data("output", "analysis-summary", "residency-patter-reanalysis.qs"))
# Summarise
residency |> 
  filter(!is.na(zone)) |>
  group_by(zone) |> 
  summarise(utils.add::basic_stats(time, na.rm = TRUE))


#### End of code. 
###########################
###########################