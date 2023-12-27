###########################
###########################
#### run-forward-sim-1.R

#### Aims
# 1) Run forward simulations

#### Prerequisites
# 1) Prepare data and parameters
# 2) Estimated storage requirements: 
#    200 GB = 700 MB * 273 individuals


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()
Sys.time()

#### Essential packages
library(dv)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(patter)
library(tictoc)

#### Load data 
dv::src()
obs_ls   <- qs::qread(here_data("input", "obs.qs"))
moorings <- readRDS(here_data("mefs", "moorings.rds"))
overlaps <- readRDS(here_data("input", "overlaps.rds"))
# bathy  <- terra::rast(here_data("spatial", "bathy.tif"))
# bset   <- terra::rast(here_data("spatial", "bset.tif"))


###########################
###########################
#### Prepare simulations 

#### Build folders (~ 2s)
tic()
sapply(sapply(obs_ls, \(d) d$folder[1]), function(folder) {
  unlink(folder, recursive = TRUE)
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
}) |> invisible()
toc()

#### Wrap SpatRasters
# TO DO
# Test speed with/without wrapping

#### (optional) Test
cl   <- 50L
test <- FALSE
if (test) {
  sel <- c("6.04-2016.dcpf", 
           "8.04-2016.acpf", 
           "13.04-2016.acdcpf")
  obs_ls <- obs_ls[sel]
  obs_ls <- lapply(obs_ls, function(obs) {
    n   <- 20L
    ind <- seq_len(n)
    if (any(c("acpf", "acdcpf") %in% c(obs$algorithm))) {
      pos_detections <- which(!sapply(obs$receiver_id, is.null))
      start <- pos_detections[1]
      start <- max(c(0, (start - ceiling(n/2))))
      ind <- start:(start + n - 1)
    }
    stopifnot(length(ind) == n)
    obs <- obs[ind, ]
    obs$timestep <- seq_len(nrow(obs))
    obs
  })
  names(obs_ls) <- obs_ls
  stopifnot(length(sel) == length(obs_ls))
  cl <- 2L
}
stopifnot(!any(is.null(sapply(obs_ls, \(elm) elm))))


###########################
###########################
#### Run simulations

#### Timings
# ETA: 11 hours if 2 hours per run and minimal parallelisation overhead
2 * length(obs_ls) / cl

gc()
tic()
pbapply::pblapply(obs_ls, cl = cl, function(obs) {

  #### Define output files
  # obs <- obs_ls[[1]]
  print(obs[1, ])
  alg        <- obs$algorithm[1]
  pff_folder <- obs$folder[1]
  log.txt    <- file.path(pff_folder, "log.txt")
  
  #### Define baseline args
  args <- list(.obs = obs,
               .bathy = terra::rast(here_data("spatial", "bathy.tif")),
               .n = 1e3, 
               .trial_kick = 1L, 
               .trial_kick_crit = 25L, 
               .trial_revert_crit = 10L, 
               .trial_revert_steps = 100L, 
               .trial_revert = 10L, 
               .save_opts = FALSE, 
               .write_opts = list(sink = pff_folder), 
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
  }
  # AC-related args
  if (alg %in% c("acpf", "acdcpf")) {
    args$.moorings           <- moorings
    args$.detection_overlaps <- overlaps
    args$.detection_kernels  <- acs_setup_detection_kernels_read()
  }
  
  #### Run algorithm 
  out_pff <- do.call(patter::pf_forward, args)
  
  #### Save outputs
  qs::qsave(out_pff, file.path(pff_folder, "out_pff.qs"))
  invisible(NULL)
  
}) |> invisible()
toc()

#### Check convergence
runs <- 
  obs_ls |> 
  rbindlist() |> 
  group_by(individual_id, block, algorithm) |> 
  slice(1L)
runs$convergence <- 
  pbapply::pbsapply(file.path(runs$folder, "out_pff.qs"), function(f) {
  qs::qread(f)$convergence
})
table(runs$convergence)


#### End of code. 
###########################
###########################