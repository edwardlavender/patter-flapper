###########################
###########################
#### prepare-patter.R

#### Aims
# 1) Prepare flapper algorithm inputs

#### Prerequisites
# 1) Obtain raw data
# 2) https://github.com/edwardlavender/flapper_appl/blob/master/R/define_global_param.R


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)

#### Essential packages
library(dv)
library(patter)
library(tictoc)

#### Load data
bathy    <- terra::rast(here_data("spatial", "bathy.tif"))
moorings <- readRDS(here_data("mefs", "moorings.rds"))

#### Local variables
overwrite <- FALSE


###########################
###########################
#### Prepare AC* inputs

#
# TO DO
# Scrap this approach and use a distances matrix?
# This is too slow, especially for large grids. 
#

#### Detection containers 
# Define containers (~50 mins)
tic()
containers <- acs_setup_detection_containers(bathy, moorings)
toc()
# Write to file for reference (~7 min)
if (overwrite) {
  pbapply::pblapply(seq_len(length(containers)), function(i) {
    container <- containers[[i]]
    if (!is.null(container)) {
      outfile <- paste0(names(containers)[i], ".tif")
      terra::writeRaster(container, here_data("input", "containers", outfile))
    }
    NULL
  }) |> invisible()
}

#### Detection overlaps
tic()
overlaps   <- acs_setup_detection_overlaps(containers, moorings)
toc()

#### Detection kernels
tic()
kernels    <- acs_setup_detection_kernels(moorings, 
                                          calc_detection_pr = acs_setup_detection_pr)
toc()

#### End of code. 
###########################
###########################