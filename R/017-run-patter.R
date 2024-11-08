###########################
###########################
#### run-patter.R

#### Aims
# 1) Run patter algorithms

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
if (!os_linux()) {
  map <- terra::rast(here_data("spatial", "bathy.tif"))
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
if (os_linux()) {
  stopifnot(!any(c("terra", "sf") %in% sort(loadedNamespaces())))
}
julia_connect()
set_seed()
set_map(here_data("spatial", "bathy.tif"))
set_vmap(.vmap = here_data("spatial", glue("vmap-{pars$pmovement$mobility[1]}.tif")))
set_model_move_components()

#### Define iterations
nrow(iteration)
iteration[, folder_xinit := file.path("data", "input", "xinit", "1.5", individual_id, month_id)]
iteration[, file_coord := file.path(folder_coord, "coord-smo.qs")]
datasets <- list(detections_by_unit = acoustics_by_unit, 
                 moorings = moorings,
                 archival_by_unit = archival_by_unit, 
                 behaviour_by_unit = behaviour_by_unit)
# Checks
all(file.exists(file.path(iteration$folder_xinit, "xinit-fwd.qs"))) |> stopifnot()
all(file.exists(file.path(iteration$folder_xinit, "xinit-bwd.qs"))) |> stopifnot()

#### Select iterations
# Start with 'best' implementations 
iteration <- iteration[sensitivity == "best", ]
# Focus on ACDC implementations (hard)
iteration <- iteration[dataset == "acdc", ]
table(table(iteration$unit_id))
# iteration <- iteration[2:.N, ]
gc()
nrow(iteration)

#### Estimate coordinates
# (TO DO) Use log.txt
lapply_estimate_coord_patter(iteration = iteration,
                             datasets = datasets, 
                             trial = FALSE, 
                             log.folder = here_data("output", "log", "real", "analysis"))

# Examine selected coords 
if (on_linux()) {
  stop("Continue on Mac (for convenience).")
}
lapply_qplot_coord(iteration, 
                   "coord-smo.qs",
                   extract_coord = function(s) s$states[sample.int(1000, size = .N, replace = TRUE), ])

#### Estimate UDs
# Time trial 
lapply_estimate_ud_spatstat(iteration = iteration[1, ], 
                            extract_coord = function(s) s$states,
                            cl = NULL, 
                            plot = FALSE)
# Implementation 
lapply_estimate_ud_spatstat(iteration = iteration, 
                            extract_coord = function(s) s$states,
                            cl = NULL, 
                            plot = FALSE)
# (optional) Examine selected UDs
lapply_qplot_ud(iteration, "spatstat", "h", "ud.tif")


#### End of code.
###########################
###########################