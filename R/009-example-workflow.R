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
# op <- options(error = recover)

#### Essential packages
library(dv)
library(data.table)
library(dtplyr)
library(dplyr)
library(ggplot2)
library(patter)
library(tictoc)

#### Load data
dv::src()
bathy     <- terra::rast(here_data("spatial", "bathy.tif"))
im        <- qs::qread(here_data("spatial", "im.qs"))
win       <- qs::qread(here_data("spatial", "win.qs"))
moorings  <- readRDS(here_data("mefs", "moorings.rds"))
acoustics <- readRDS(here_data("mefs", "acoustics.rds"))
archival  <- readRDS(here_data("mefs", "archival.rds"))
overlaps  <- readRDS(here_data("input", "overlaps.rds"))
kernels   <- acs_setup_detection_kernels_read()

#### Local pars
seed <- 1L
run    <- FALSE
manual <- run


###########################
###########################
#### Set up analysis

#### Collate observations
# Limits
xlim <- as.POSIXct(c("2016-07-01 00:00:00", "2016-08-01 23:58:00"), tz = "UTC")
# Acoustics
acc <- 
  acoustics |>
  filter(individual_id == 25) |> 
  filter(timestamp >= xlim[1]) |> 
  filter(timestamp <= xlim[2]) |> 
  as.data.table()
# Archival data (if applicable)
arc <- 
  archival |>
  filter(individual_id == 25) |> 
  filter(timestamp >= xlim[1]) |> 
  filter(timestamp <= xlim[2]) |> 
  arrange(timestamp) |>
  # Define 'resting': 0L = resting; 1L = not resting;
  mutate(va = abs(serial_difference(depth)), 
         state = if_else(va <= 0.25, true = 0L, false = 1L)) |>
  as.data.table()
arc$va <- NULL
arc$state[nrow(arc)] <- 1L
# Collate observations 
obs <- acs_setup_obs(acc, arc, 
                     .trim = FALSE,
                     .step = "2 mins", 
                     .mobility = 500, 
                     .detection_range = 750)
obs[, state := arc$state[match(timestamp, arc$timestamp)]]
# Define behavioural state
# Check obs
obs
range(obs$timestamp)
stopifnot(any(obs$detection == 1L))
stopifnot(!any(is.na(obs$depth)))
# Visualise obs
plot(arc$timestamp, arc$depth*-1, 
     xlim = xlim, ylim = c(-225, 0),
     type = "n")
s <- seq(nrow(arc) - 1)
arrows(arc$timestamp[s], arc$depth[s] * -1, 
       arc$timestamp[s + 1L], arc$depth[s + 1] * -1, 
       col = factor(arc$state), length = 0, lwd = 0.5)
points(acc$timestamp, rep(0L, nrow(acc)), col = "red")

#### Define algorithm 
algs <- c("acpf", "acdcpf", "dcpf")
alg <- algs[2]


###########################
###########################
#### Forward filter


###########################
#### Set up 

#### Define directories
pff_folder <- here_data("example", "forward", alg, "output")
log.txt    <- here_data("example", "forward", alg, "log.txt")
# Optionally wipe old directories
if (FALSE) {
  unlink(log.txt)
  unlink(pff_folder, recursive = TRUE)
}
# (Re)build
if (!dir.exists(pff_folder)) {
  dir.create(pff_folder, recursive = TRUE)
}

#### Define baseline args
args <- list(.obs = obs,
             .bathy = terra::rast(here_data("spatial", "bathy.tif")),
             .n = 1e3, 
             .trial_kick = 1L, 
             .trial_kick_crit = 25L, 
             .trial_revert_crit = 10L, 
             .trial_revert_steps = 100L, 
             .trial_revert = 10L,
             .record_opts = list(save = FALSE,
                                 sink = pff_folder, 
                                 cols = c("timestep", "cell_past", "cell_now", "x_now", "y_now")),
             .progress = FALSE,
             .verbose = TRUE, 
             .txt = log.txt) 

#### Define algorithm-specific args 
# DC-related args 
if (alg %in% c("dcpf", "acdcpf")) {
  # Define origin starting point 
  args$.origin    <- dc_origin(.bathy = args$.bathy, 
                               .depth = obs$depth[1], 
                               .calc_depth_error = calc_depth_error)
  # Define DC model
  bset <- terra::rast(here_data("spatial", "bset.tif"))
  args$.update_ac <- function(.particles, .bathy, .obs, .t, 
                              .bset = bset, .calc_depth_error = calc_depth_error) {
    update_ac(.particles = .particles, .bathy = .bathy, .obs = .obs, .t = .t, 
              .bset = .bset, .calc_depth_error = .calc_depth_error)
  }
  # Define behaviourally dependent movement model (via .rpropose and .dpropose)
  # * TO DO
}
# AC-related args
if (alg %in% c("acpf", "acdcpf")) {
  args$.moorings           <- moorings
  args$.detection_overlaps <- overlaps
  args$.detection_kernels  <- acs_setup_detection_kernels_read()
}
if (alg == "acpf") {
  # Define .rpropose and .dpropose
  # * Use default movement models
}

#### (optional) TO DO
# Downgrade movement jumps using shortest distances model


###########################
#### Forward run 

#### Timings
# * ACPF: ~107 mins
# * DCPF: 
# * ACDCPF: 

if (TRUE) {
  tic()
  set.seed(seed)
  out_pff <- do.call(patter::pf_forward, args)
  toc()
  # beepr::beep(10L)
}


###########################
#### Outputs

#### Output file size (MB)
hs <- file.size(list.files(file.path(pff_folder, "history"), full.names = TRUE))
ds <- file.size(list.files(file.path(pff_folder, "diagnostics"), full.names = TRUE))
size <- sum(c(hs, ds))
size/1e6 # MB

#### Output diagnostics 
if (manual) {
  diag <- pf_forward_diagnostics(pff_folder)
  tail(diag, 100)
}


###########################
#### Validation 

#### Validate the number of particles at each time step (~2 s)
pff_folder_h <- file.path(pff_folder, "history")
tic()
np <- 
  pff_folder_h |> 
  arrow::open_dataset() |> 
  group_by(timestep) |> 
  dplyr::count() |> 
  arrange(timestep) |>
  collect()
which(np$n != 1000L)
stopifnot(all(np$n == 1e3L))
toc()

#### TO DO
# 


###########################
###########################
#### Backward sampler

#### Define directories
pfb_folder  <- here_data("example", "backward", alg)
pfbk_folder <- file.path(pfb_folder, "killer", "output")
log.txt     <- file.path(pfb_folder, "killer", "log.txt")
if (FALSE) {
  unlink(pfbk_folder, recursive = TRUE)
  unlink(log.txt)
}
if (!dir.exists(pfbk_folder)) {
  dir.create(pfbk_folder, recursive = TRUE)
}

#### Implement backward killer 
# ACPF: ~140 s
if (run) {
  tic()
  out_pfbk <- pf_backward_killer(pf_setup_files(pff_folder_h), 
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
}


###########################
###########################
#### Backward sampler

#### Define directories 
pfbs_folder <- file.path(pfb_folder, "sampler", "output")
log.txt     <- file.path(pfb_folder, "sampler", "log.txt")

#### Run sampler 
files <- pf_setup_files(pff_folder_h)
files <- files[(length(files) - 10L):length(files)]
if (run) {
  tic()
  set.seed(seed)
  out_pfbs <- pf_backward_sampler(files, 
                                  .step_dens = dstep, lonlat = FALSE,
                                  .write_history = list(sink = pfbs_folder), 
                                  .txt = log.txt)
  toc()
}


###########################
###########################
#### Maps

#### Define input folder
folder <- pfbk_folder
# folder <- pfbs_folder

#### Collect coordinates (~5 s)
# * We collect coordinates and use pf_map_dens()
# * This is faster than pf_map_pou() plus pf_map_dens()
tic()
pxy <- pf_coords(folder, .bathy = bathy)
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
                   here_data("example", "map", "acpf", "dens.tif"), 
                   overwrite = TRUE)
toc()


#### End of code
###########################
###########################