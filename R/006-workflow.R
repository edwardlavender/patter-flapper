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
dv::clear()
options(error = recover)

#### Essential packages
library(dv)
library(data.table)
library(dtplyr)
library(dplyr)
library(ggplot2)
library(patter)
library(tictoc)

#### Load data
source(here_r("002-define-helpers.R"))
bathy     <- terra::rast(here_data("spatial", "bathy.tif"))
im        <- qs::qread(here_data("spatial", "im.qs"))
win       <- qs::qread(here_data("spatial", "win.qs"))
moorings  <- readRDS(here_data("mefs", "moorings.rds"))
acoustics <- readRDS(here_data("mefs", "acoustics.rds"))
archival  <- readRDS(here_data("mefs", "archival.rds"))
overlaps  <- readRDS(here_data("input", "overlaps.rds"))
kernels   <- acs_setup_detection_kernels_read()

#### Local pars
set.seed(1)
run    <- FALSE
manual <- run


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
#### Forward filter

#### TO DO
# Write behavioural switching model
# Downgrade movement jumps using shortest distances model

#### Define directories
log.txt    <- here_data("example", "acpf", "forward", "log.txt")
pff_folder <- here_data("example", "acpf", "forward", "output")
if (FALSE) {
  unlink(log.txt)
  unlink(pff_folder, recursive = TRUE)
  dir.create(pff_folder, recursive = TRUE)
}

#### Forward run 
# ACPF: ~192 mins
tic()
if (run) {
  out_pff <- pf_forward(obs,
                        .bathy = bathy,
                        .moorings = moorings,
                        .detection_overlaps = overlaps,
                        .detection_kernels = kernels,
                        .n = 1e3L,
                        # .update_ac = update_ac,
                        .trial_kick = 1L,
                        # Use .trial_kick_crit to control directed sampling threshold 
                        .trial_kick_crit = 25L, 
                        .trial_revert_crit = 10L, 
                        .trial_revert_steps = 100L, 
                        .trial_revert = 10L,
                        .save_opts = TRUE,
                        .write_opts = list(sink = pff_folder), 
                        .verbose = TRUE, .txt = log.txt)
  toc()
  # beepr::beep(10L)
}

#### Diagnostics
if (manual) {
  diag <- pf_forward_diagnostics(pff_folder)
  tail(diag, 100)
}


###########################
###########################
#### Backward killer

#### Define directories
pfb_folder  <- here_data("example", "acpf", "backward")
pfbk_folder <- file.path(pfb_folder, "killer", "output")
log.txt     <- here_data("example", "acpf", "backward", "killer", "log.txt")

#### Implement pruning 
# ACPF: ~140 s
if (run) {
  tic()
  out_pfbk <- pf_backward_killer(pf_setup_files(file.path(pff_folder, "history")), 
                                 .write_history = list(sink = pfbk_folder), 
                                 .txt = log.txt)
  toc()
}

#### Diagnostics 
if (manual) {
  # Extract diagnostics 
  # * 1 CPU: ~38 s
  # * 10 CPU (fork): 7 s (with chunks)
  pfbk_diag <- pf_backward_killer_diagnostics(pfbk_folder)
  # Examine diagnostic traces
  plot(pfbk_diag$timestep, pfbk_diag$n_u, type = "l")
  plot(pfbk_diag$timestep, pfbk_diag$n_u, ylim  = c(0, 100), type = "l")
  plot(pfbk_diag$timestep, pfbk_diag$ess, type = "l")
  # TO DO
  # * Resolve issue in pf_forward() (?)
  # * How can we have more unique cells than particles?
  pfbk_diag[n_u > 1000, ]
  arrow::read_parquet(file.path(pfbk_folder, "3252.parquet"))
}


###########################
###########################
#### Backward sampler

#### Define directories 
pfbs_folder <- file.path(pfb_folder, "sampler", "output")
log.txt     <- here_data("example", "acpf", "backward", "killer", "log.txt")

#### Prepare sampler
# Identify distinct cells (~3.2 mins, 1 CPU)
tic()
pfbd_folder <- file.path(pfb_folder, "sampler", "distinct")
pf_distinct(.history = file.path(pff_folder, "history"), 
            .write_opts = list(sink = pfbd_folder))
toc()

#### Run sampler 
if (run) {
  tic()
  out_pfbs <- pf_backward_sampler(pf_setup_files(pff_folder), 
                                  .step_dens = dstep, 
                                  .write_history = list(sink = pfbs_folder), 
                                  .txt = log.txt)
  toc()
}


###########################
###########################
#### Maps

#### TO DO
# * Update to handle multiple intputs
# * pfbk_folder, pfbs_folder

#### Collect coordinates (~5 s)
# * We collect coordinates and use pf_map_dens()
# * This is faster than pf_map_pou() plus pf_map_dens()
tic()
pxy <- pf_coords(pfbk_folder, .bathy = bathy)
toc()
if (manual) {
  # Number of unique cells sampled
  length(unique(pxy$cell_id))
  # Unique samples through time
  pxy |> 
    group_by(timestep) |> 
    summarise(n = n()) |> 
    ggplot() + 
    geom_line(aes(timestep, n))
}

#### Density (~70 s)
tic()
dens <- pf_map_dens(.xpf =  bathy, 
                    .im = im, 
                    .owin = win, 
                    .coord = pxy
                    )
if (manual) {
  # Add unique particle samples
  pxyu <- pxy[!duplicated(cell_id), ]
  points(pxyu$cell_x, pxyu$cell_y, 
         pch = ".", 
         col = scales::alpha("darkgrey", 0.7))
}
toc()

#### Write maps to file (~11 s)
tic()
terra::writeRaster(dens, 
                   here_data("example", "acpf", "map", "dens.tif"), 
                   overwrite = TRUE)
toc()


#### End of code
###########################
###########################