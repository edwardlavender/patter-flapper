###########################
###########################
#### run-rsp.R

#### Aims
# 1) Run RSP algorithm 

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
library(patter)
dv::src()

#### Load data 
iteration <- qs::qread(here_data("input", "iteration", "rsp.qs"))
moorings  <- qs::qread(here_data("input", "mefs", "moorings.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))


###########################
###########################
#### Run algorithm 

#### (optional) Testing
test <- TRUE
if (test) {
  iteration <- iteration[1:2]
} 

#### Estimate coordinates (RSPs)
nrow(iteration)
moorings[, receiver_gamma := 1500]
datasets <- list(detections_by_unit = acoustics_by_unit, moorings = moorings)
lapply_estimate_coord_rsp(iteration = iteration, datasets = datasets)

#### Estimate DBBMMs
iteration[, file_coord := file.path(folder_coord, "coord.qs")]
lapply_estimate_ud_dbbmm(iteration = iteration, 
                         cl = NULL, 
                         plot = FALSE)


#### End of code.
###########################
###########################