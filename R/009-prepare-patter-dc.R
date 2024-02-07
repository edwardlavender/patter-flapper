###########################
###########################
#### prepare-patter-dc.R

#### Aims
# 1) Prepare *DCPF algorithm components:
#    * Depth-error windows

#### Prerequisites
# 1) Obtain raw data
# 2) https://github.com/edwardlavender/flapper_appl/blob/master/R/define_global_param.R


###########################
###########################
#### Set up 

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(dv)
library(tictoc)

#### Load data
dv::src()
bathy <- terra::rast(here_data("spatial", "bathy.tif"))
bset  <- terra::rast(here_data("spatial", "bset.tif"))


###########################
###########################
#### Define depth-error window

# ~6.5 mins
if (FALSE) {
  tic()
  ewin <- dc_ewindow(.bathy = bathy, .bset = bset)
  writeRasterLs(x = ewin, 
                folder = here_data("input", "depth-window"), 
                index = FALSE)
  toc()
}


#### End of code. 
###########################
###########################