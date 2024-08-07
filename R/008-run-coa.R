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
library(dv)
library(patter)
library(lubridate)
dv::src()

#### Load data 
iteration <- qs::qread(here_data("input", "iteration", "coa.qs"))
moorings  <- qs::qread(here_data("input", "mefs", "moorings.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))


###########################
###########################
#### Run algorithm 

#### (optional) Testing
test <- FALSE
if (test) {
  iteration <- iteration[1:2]
} 

#### Estimate coordinates (COAs): ~19 s
nrow(iteration)
datasets <- list(detections_by_unit = acoustics_by_unit, moorings = moorings)
lapply_estimate_coord_coa(iteration = iteration, datasets = datasets)

#### Time trials (to estimate 1 density surface):
# Test sigma = NULL with npixel (modify src/options.R):
# * 50 pixel: 2 s, too coarse
# * 250 pixel: 42 s, too coarse
# * 400 pixel: 106 s, a bit too coarse
# * 500 pixel, 165 s, adequate 
# > See fig/spatstat for images
if (FALSE) {
  lapply_estimate_ud_spatstat(iteration = iteration[1, ], cl = NULL)
}
# Estimated total duration for all iterations on 1 cl (hours):
nrow(iteration) * 165 / 60 / 60

#### Estimate UDs
iteration[, file_coord := file.path(folder_coord, "coord.qs")]
lapply_estimate_ud_spatstat(iteration = iteration, 
                            extract_coord = NULL,
                            cl = NULL, 
                            plot = FALSE)


#### End of code.
###########################
###########################