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
library(data.table)
library(dtplyr)
library(dplyr)
library(patter)
library(tictoc)

#### Load data
source(here_r("002-define-helpers.R"))
bathy     <- terra::rast(here_data("spatial", "bathy.tif"))
moorings  <- readRDS(here_data("mefs", "moorings.rds"))
acoustics <- readRDS(here_data("mefs", "acoustics.rds"))
archival  <- readRDS(here_data("mefs", "archival.rds"))


###########################
###########################
#### Select time series for analysis 

#### Collate observations
# Acoustics
acc <- 
  acoustics |>
  filter(individual_id == 25) |> 
  filter(timestamp >= as.POSIXct("2016-07-01")) |> 
  filter(timestamp <= as.POSIXct("2016-08-01")) |> 
  as.data.table()
# Archival data (if applicable)
arc <- 
  archival |>
  filter(individual_id == 25) |> 
  filter(timestamp >= as.POSIXct("2016-07-01")) |> 
  filter(timestamp <= as.POSIXct("2016-08-01")) |> 
  as.data.table()
# Collate observations 
obs <- acs_setup_obs(acc, arc, 
                     .step = "2 mins", 
                     .mobility = 500, 
                     .detection_range = 750)
# Add depth errors (if applicable)
depth_error <- calc_depth_error(obs$depth)
obs[, depth_shallow := depth + depth_error[1, ]]
obs[, depth_deep := depth + depth_error[2, ]]


###########################
###########################
#### Implement forward filter--backward sampler

#### TO DO
# Write behavioural switching model
# Downgrade movement jumps using shortest distances model

#### Forward run 
tic()
log.txt <- here_data("example", "forward", "log.txt")
pff_folder <- here_data("example", "forward", "output")
out_pff <- pf_forward_2(obs, 
                        .bathy = bathy, 
                        .moorings = moorings, 
                        .detection_overlaps = overlaps, 
                        .detection_kernels = kernels, 
                        .update_ac = update_ac, 
                        .kick = pf_kick, 
                        .n = 1e3, 
                        .save_history = FALSE, 
                        .write_history = list(sink = pff_folder), 
                        .verbose = TRUE, .txt = log.txt)
toc()

#### Prepare backward sampler
# (optional) TO DO
# * Trial in memory computations
tic()
log.txt    <- here_data("example", "backward", "input", "log.txt")
pfd_folder <- here_data("example", "backward", "input")
pf_backward_dens(pff_folder, 
                 .step_dens = step_dens, 
                 .in_memory = FALSE, 
                 .store = pfd_folder,
                 .verbose = TRUE, .txt = log.txt)
toc()

#### Run backwards sample
tic()
log.txt    <- here_data("example", "backward", "log.txt")
pfd_folder <- here_data("example", "backward", "input", "density")
pfb_folder <- here_data("example", "backward", "output")
input      <- pf_setup_files(pff_folder)
pf_backward_p(input, 
              .step_dens = step_dens_read, 
              .density = pfd_folder,
              .save_history = FALSE, 
              .write_history = list(sink = pfb_folder), 
              .verbose = TRUE, .txt = log.txt)
toc()


###########################
###########################
#### Analyse outputs

# TO DO
# * Use a low resolution grid for weights in estimation


#### End of code
###########################
###########################