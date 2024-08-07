###########################
###########################
#### run-coa.R

#### Aims
# 1) Run COA-KUD algorithm 

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
bathy     <- terra::rast(here_data("spatial", "bathy.tif"))
iteration <- qs::qread(here_data("input", "iteration", "coa.qs"))
moorings  <- qs::qread(here_data("input", "mefs", "moorings.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))


###########################
###########################
#### Run algorithm 

#### Set up
nrow(iteration)
iteration[, file_coord := file.path(folder_coord, "coord.qs")]
datasets <- list(detections_by_unit = acoustics_by_unit, moorings = moorings)

#### (optional) Testing
test <- TRUE
if (test) {
  iteration <- iteration[1:2]
} 

#### Estimate coordinates
# Time trial
lapply_estimate_coord_coa(iteration = iteration[1, ], datasets = datasets)
# Implementation (~19 s)
lapply_estimate_coord_coa(iteration = iteration, datasets = datasets)
# (optional) Examine selected coords
lapply_qplot_coord(iteration, "coord.qs")

#### Estimate UDs
# Time trial 
lapply_estimate_ud_spatstat(iteration = iteration[1, ], 
                            extract_coord = NULL,
                            cl = NULL, 
                            plot = FALSE)
# Implementation 
lapply_estimate_ud_spatstat(iteration = iteration, 
                            extract_coord = NULL,
                            cl = NULL, 
                            plot = FALSE)
# (optional) Examine selected UDs
lapply_qplot_ud(iteration, "spatstat", "h", "ud.tif")


#### End of code.
###########################
###########################