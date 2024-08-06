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
library(dv)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)
library(JuliaCall)
library(patter)
dv::src()

#### Load data 
map       <- terra::rast(here_data("spatial", "bathy.tif"))
iteration <- qs::qread(here_data("input", "iteration", "patter.qs"))
moorings  <- qs::qread(here_data("input", "mefs", "moorings.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))
archival_by_unit  <- qs::qread(here_data("input", "archival_by_unit.qs"))


###########################
###########################
#### Run algorithm 

#### Connect to Julia 
julia_connect()
set_seed()
set_map(map)
JuliaCall::julia_command(ModelObsAcousticContainer)
JuliaCall::julia_command(ModelObsAcousticContainer.logpdf_obs)

#### (optional) Testing
test <- TRUE
if (test) {
  iteration <- iteration[1:2L, ]
} 

#### Estimate coordinates (patter)
nrow(iteration)
iteration[, index := 1:.N]
datasets <- list(detections_by_unit = acoustics_by_unit, 
                 moorings = moorings,
                 archival_by_unit = archival_by_unit)
lapply_estimate_coord_patter(iteration = iteration, 
                             datasets = datasets)

#### Estimate UDs
iteration[, file_coord := file.path(folder, "coord-smo.qs")]
lapply_estimate_ud_spatstat(iteration = iteration, 
                            extract_coord = function(x) x$states,
                            cl = NULL, 
                            plot = FALSE)


#### End of code.
###########################
###########################