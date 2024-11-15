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
if (!patter:::os_linux()) {
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
if (patter:::os_linux()) {
  stopifnot(!any(c("terra", "sf") %in% sort(loadedNamespaces())))
}
julia_connect()
set_seed()
set_map(here_data("spatial", "bathy.tif"))
set_model_move_components()

#### Clean up
if (FALSE) {
  unlink(iteration$folder_coord, recursive = TRUE)
  unlink(iteration$folder_ud, recursive = TRUE)
  list.files(iteration$folder_coord, recursive = TRUE)
  list.files(iteration$folder_ud, recursive = TRUE)
  dirs.create(iteration$folder_coord)
  dirs.create(iteration$folder_ud)
}

#### Prepare iterations/datasets
iteration[, folder_xinit := file.path("data", "input", "xinit", "1.5", individual_id, month_id)]
iteration[, file_coord := file.path(folder_coord, "coord-smo.qs")]
datasets <- list(detections_by_unit = acoustics_by_unit, 
                 moorings           = moorings,
                 archival_by_unit   = archival_by_unit, 
                 behaviour_by_unit  = behaviour_by_unit)

#### Review inputs
# Check iterations
# * Each unit_id is a individual/month combination
# * For each unit_id, we have AC/DC/ACDC algorithm implementations
# * For each algorithm implementation, we have different parameterisations (sensitivity)
head(iteration)
# Check datasets 
# * unit_ids were defined as _all possible_ individual-month combinations
# * (see process-data-mefs.R)
# -> We have 154 possible unit ids (individual-month) combinations
ids  <- unique(iteration$individual_id)
mons <- mmyy(seq(as.Date("2016-04-01"), as.Date("2017-05-31"), by = "months"))
nrow(CJ(ids, mons))
# * We have 48 'realised' unit_ids 
# -> These are the individual/months for which we have observations
sort(unique(iteration$unit_id))
length(unique(iteration$unit_id))
# * We have 154 elements in each dataset (one for each _possible_ unit_id)
length(datasets$detections_by_unit)
length(datasets$archival_by_unit)  
length(datasets$behaviour_by_unit)
# * For all AC/ACDC unit_id, we need detections dataset
sapply(unique(iteration$unit_id[iteration$dataset %in% c("ac", "acdc")]), function(id) {
  !is.null(datasets$detections_by_unit[[id]])
}) |> all() |> stopifnot()
# * For all DC/ACDC unit_id, we need archival dataset
sapply(unique(iteration$unit_id[iteration$dataset %in% c("dc", "acdc")]), function(id) {
  !is.null(datasets$archival_by_unit[[id]])
}) |> all() |> stopifnot()
# * For all AC/DC/ACDC unit_id, we need behaviour dataset
sapply(unique(iteration$unit_id), function(id) {
  !is.null(datasets$behaviour_by_unit[[id]])
}) |> all() |> stopifnot()
# Check folder_xinit
all(file.exists(file.path(iteration$folder_xinit, "xinit-fwd.qs"))) |> stopifnot()
all(file.exists(file.path(iteration$folder_xinit, "xinit-bwd.qs"))) |> stopifnot()

#### Select batch
# As for the simulations, we batch by mobility
batch     <- pars$pmovement$mobility[1]
iteration <- iteration[mobility == batch, ]
set_vmap(.vmap = here_data("spatial", glue("vmap-{batch}.tif")))

#### Select iterations
# Select iterations 
# * Select one AC/DC/ACDC run for testing 
nrow(iteration)
# iteration <- 
#   iteration |> 
#   filter(sensitivity == "best") |>
#   group_by(dataset) |> 
#   slice(1L) |>
#   as.data.table()
# * Select by dataset/sensitivity
iteration <- iteration[dataset == "acdc", ]
iteration <- iteration[sensitivity == "best", ]

#### Estimate coordinates
# * Convergence trials complete (see /log/real/trials/log-summary.txt)
# * Batch 1, AC (all): SIA-USER024-P
# * Batch 1, DC (all): MCC02XT0AZJGH5
# * Batch 2, ACDC    : TO DO
gc()
nrow(iteration)
lapply_estimate_coord_patter(iteration  = iteration,
                             datasets   = datasets, 
                             trial      = FALSE, 
                             log.folder = here_data("output", "log", "real"), 
                             log.txt    = NULL)

#### Estimate UDs
if (patter:::os_linux()) {
  
  # Examine selected coord
  # lapply_qplot_coord(iteration, 
  #                    "coord-smo.qs",
  #                    extract_coord = function(s) s$states[sample.int(1000, size = .N, replace = TRUE), ])
  
  #### Estimate UDs
  # Time trial 
  lapply_estimate_ud_spatstat(iteration     = iteration[1, ], 
                              extract_coord = function(s) s$states,
                              cl            = NULL, 
                              plot          = FALSE)
  # Implementation 
  lapply_estimate_ud_spatstat(iteration     = iteration, 
                              extract_coord = function(s) s$states,
                              cl            = NULL, 
                              plot          = FALSE)
  # (optional) Examine selected UDs
  lapply_qplot_ud(iteration, "spatstat", "h", "ud.tif")
}


#### End of code.
###########################
###########################