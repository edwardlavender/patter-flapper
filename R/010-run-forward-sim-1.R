###########################
###########################
#### run-forward-sim-1.R

#### Aims
# 1) Run forward simulations

#### Prerequisites
# 1) Prepare data and parameters


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
library(ggplot2)
library(patter)
library(tictoc)

#### Load data 
src()
skateids  <- readRDS(here_data("mefs", "skateids.rds"))
acoustics <- readRDS(here_data("mefs", "acoustics.rds"))
archival  <- readRDS(here_data("mefs", "archival.rds"))
bathy     <- terra::rast(here_data("spatial", "bathy.tif"))


###########################
###########################
#### Run simulations 

#### Prepare directories
# data/output/forward/{individual}/{mmyy}/{algorithm}/output/
# > log.txt
# > history/
# > diagnostics 

#### Run simulations
tic()
pbapply::pblapply(obs, cl = NULL, function(o) {
  
  
  
})
toc()


#### End of code. 
###########################
###########################