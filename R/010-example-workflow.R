###########################
###########################
#### example-workflow.R

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
op <- options(error = recover)

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
bset      <- terra::rast(here_data("spatial", "bset.tif"))
im        <- qs::qread(here_data("spatial", "im.qs"))
win       <- qs::qread(here_data("spatial", "win.qs"))
moorings  <- readRDS(here_data("mefs", "moorings.rds"))
acoustics <- readRDS(here_data("mefs", "acoustics.rds"))
archival  <- readRDS(here_data("mefs", "archival.rds"))
overlaps  <- readRDS(here_data("input", "overlaps.rds"))
kernels   <- acs_setup_detection_kernels_read()
ewin      <- readRasterLs(here_data("input", "depth-window"))

#### Local pars
seed <- 1L
run    <- FALSE
manual <- run


###########################
###########################
#### Set up analysis

#### Collate observations
# Limits
xlim <- as.POSIXct(c("2016-07-01 00:00:00", "2016-08-01 23:58:00"), 
                   tz = "UTC")
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
# Build data list
# * We add a .$spatial$origin element optionally below (for DCPF/ACDCPF algorithms)
dlist <- pat_setup_data(.acoustics = acc, 
                        .moorings = moorings, 
                        .archival = arc, 
                        .bathy = bathy, 
                        .lonlat = FALSE)
# Include AC likelihood terms (required for ACPF, ACDCPF)
dlist$algorithm$detection_overlaps <- overlaps
dlist$algorithm$detection_kernels  <- kernels
# Include DC likelihood terms (required for DCPF, ACDCPF)
dlist$spatial$bset      <- bset
dlist$algorithm$ewindow <- ewin
# Define observations timeline
obs <- pf_setup_obs(.dlist = dlist,
                    .trim = FALSE,
                    .step = "2 mins", 
                    .mobility = 500, 
                    .receiver_range = 750)
# Define behavioural states
obs[, state := arc$state[match(timestamp, arc$timestamp)]]
# Define obs columns & obs-related dlist elements (for acs_setup_container_acdc())
dlist$algorithm$pos_detections <- which(!sapply(obs$receiver_id, is.null))
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

#### Define origin
# This is included in dlist below, if necessary
origin <- dc_origin(.ewindow = dlist$algorithm$ewindow, .depth = obs$depth[1])

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
if (run) {
  unlink(log.txt)
  unlink(pff_folder, recursive = TRUE)
}
# (Re)build
if (!dir.exists(pff_folder)) {
  dir.create(pff_folder, recursive = TRUE)
}

#### Define data list
if (alg %in% c("dcpf", "acdcpf")) {
  dlist$spatial$.origin <- origin
} else {
  dlist$spatial$.origin <- NULL
}

#### Define args
# Define baseline forward and backward arguments 
record <- 
  pf_opt_record(
    .save = FALSE,
    .sink = pff_folder, 
    .cols = c("timestep", "cell_past", "cell_now", "x_now", "y_now", "lik")
  )
# Define baseline forward arguments
# * TO DO
# * Define behaviourally dependent movement model (via .rpropose and .dpropose)
args <- list(.obs = obs,
             .dlist = dlist,
             .n = 1e3, 
             .trial = 
               pf_opt_trial(
                 # Kick once 
                 .trial_kick = 1L, 
                 # Initiate directed sampling when there are < 25 cells
                 .trial_sampler = 100L, 
                 # Revert when there are < 5 grid cells, by 50 steps, 10 times
                 .trial_revert_crit = 5L, 
                 .trial_revert_steps = 50L, 
                 .trial_revert = 10L
               ),
             .record = record,
             .control = pf_opt_control(.sampler_batch_size = 1000L),
             .verbose = log.txt) 
# Define algorithm-specific likelihoods
if (alg == "acpf") {
  args$.likelihood <- list(pf_lik_ac = pf_lik_ac, 
                           acs_filter_land = acs_filter_land, 
                           acs_filter_container = acs_filter_container)
} else if (alg == "dcpf") {
  args$.likelihood <- list(pf_lik_dc = pf_lik_dc_2)
} else if (alg == "acdcpf") {
  args$.likelihood <- list(pf_lik_dc = pf_lik_dc_2, 
                           pf_lik_ac = pf_lik_ac, 
                           acs_filter_container = acs_filter_container)
}

#### (optional) TO DO
# Downgrade movement jumps using shortest distances model


###########################
#### Forward run 

#### Timings
# * ACPF: ~107 mins
# * DCPF: 
# * ACDCPF: 

if (run) {
  tic()
  ss()
  out_pff <- do.call(patter::pf_forward, args)
  toc()
  # beepr::beep(10L)
}

# To debug convergence issues, see ./R/supporting/convergence/.

###########################
#### Outputs

#### Output file size (MB)
hs <- pf_files_size(pff_folder, .folder = "history")
ds <- pf_files_size(pff_folder, .folder = "diagnostics")
sum(c(hs, ds))

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
record$.sink <- pfbk_folder

#### Implement backward killer 
# ACPF: ~140 s
if (run) {
  tic()
  out_pfbk <- pf_backward_killer(.history = pf_files(pff_folder_h), 
                                 .record = record,
                                 .verbose = log.txt)
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
files <- pf_files(pff_folder_h)
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