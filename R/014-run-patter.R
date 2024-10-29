###########################
###########################
#### run-patter.R

#### Aims
# 1) Run patter algorithms

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
dv::src()

#### Load data 
pars              <- qs::qread(here_data("input", "pars.qs"))
map               <- terra::rast(here_data("spatial", "bathy.tif"))
vmap              <- terra::rast(here_data("spatial", glue("vmap-{pars$pmovement$mobility[1]}.tif")))
iteration         <- qs::qread(here_data("input", "iteration", "patter.qs"))
moorings          <- qs::qread(here_data("input", "mefs", "moorings.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))
archival_by_unit  <- qs::qread(here_data("input", "archival_by_unit.qs"))
behaviour_by_unit <- qs::qread(here_data("input", "behaviour_by_unit.qs"))


###########################
###########################
#### Run algorithm 

#### Julia Set up
julia_connect()
set_seed()
set_map(map)
set_vmap(.vmap = vmap)
julia_command(ModelMoveFlapper)

#### Define iterations
nrow(iteration)
iteration[, file_coord := file.path(folder_coord, "coord-smo.qs")]
datasets <- list(detections_by_unit = acoustics_by_unit, 
                 moorings = moorings,
                 archival_by_unit = archival_by_unit, 
                 behaviour_by_unit = behaviour_by_unit)

#### (optional) Testing
test <- TRUE
if (test) {
  
  # Visualise time series to pick an example individual with good time series
  if (FALSE) {
    nid <- length(unique(iteration$unit_id))
    png(here_fig("time-series-acdc.png"),
        height = 10, width = 15, units = "in", res = 600)
    pp <- par(mfrow = par_mf(nid), oma = c(2, 2, 2, 2), mar = c(1, 1, 1, 1))
    cl_lapply(unique(iteration$unit_id), function(id) {
      # id <- 2
      acc <- datasets$detections_by_unit[[id]]
      arc <- datasets$archival_by_unit[[id]]
      plot(arc$timestamp, arc$depth * -1, ylim = c(-250, 0), 
           xlab = "", ylab = "",
           type = "l", main = id)
      points(acc$timestamp, rep(0, nrow(acc)), col = "red")
    })
    par(pp)
    dev.off()
  }
  
  # Select AC/DC/ACDC implementations for an example individual with good time series
  iteration <- iteration[unit_id == 119 & sensitivity == "best", ] # fails around ~560
  # iteration <- iteration[unit_id == 55 & sensitivity == "best", ]  # fails around ~15952s, works with 100,000 particles
  iteration <- iteration[dataset == "acdc", ]
  
  # Restrict np for speed 
  # * This must be >= n_record which is 1000L
  iteration[, np := 50000L]
  
  # (optional) Trial adjusted parameters
  # * Large sigma, large depth_deep_eps -> weak influence of depths
  # * Large sigma, small depth_deep_eps -> uniform model with truncation
  # * Small sigma, large depth_deep_eps -> remove truncation
  iteration[, depth_sigma := 10000]
  iteration[, depth_deep_eps := 20]
  curve(dnorm(x, 0, 20), from = 0, to = 350)
  curve(dnorm(x, 200, 20), from = 0, to = 350)
  
  # (optional) Turn off smoothing for convergence trials
  iteration[, smooth := FALSE, ]
} 

#### Select iterations
# Start with 'best' implementations 
iteration <- iteration[sensitivity == "best", ]
# Focus on ACDC implementations (hard)
iteration <- iteration[dataset == "acdc", ]
table(table(iteration$unit_id))
# iteration <- iteration[2:.N, ]

#### Estimate coordinates
gc()
nrow(iteration)
rm(map, vmap)
# Set up log.txt file (on unix only)
if (FALSE & on_unix()) {
  dirs.create(here_data("output", "log", "analysis"))
  log.txt <- here_data("output", "log", "analysis", "patter-log.txt")
  log.txt <- file(log.txt, open = "wt")
  sink(log.txt)
  sink(log.txt, type = "message")
}
# Estimate coordinates 
Sys.time()
lapply_estimate_coord_patter(iteration = iteration, datasets = datasets)
Sys.time()
# Close sink
if (on_unix()) {
  sink()
  sink(type = "message") 
}
# Examine selected coords 
lapply_qplot_coord(iteration, 
                   "coord-smo.qs",
                   extract_coord = function(s) s$states[sample.int(1000, size = .N, replace = TRUE), ])

#### Record (by index)
# * Currently on row 86...
iteration[1:20, .(index, dataset, sensitivity, np)]
# Forward failure   : 11, 23 (DC), 28, 45, 62, 
# > ACDC unless labelled otherwise
# Backward failure  : 6 (DC), 8:10 (DCs), 57 (DC), 
# > ACDC unless labelled otherwise
# OK                : 1:5, 7, 18, 35, 40 (DC), 52, 69, 74 (DC)
# > AC unless labelled otherwise 

# * Multiple trials do help convergence, in some instances
# e.g., 74

#### Estimate UDs
# Time trial 
lapply_estimate_ud_spatstat(iteration = iteration[1, ], 
                            extract_coord = function(s) s$states,
                            cl = NULL, 
                            plot = FALSE)
# Implementation 
lapply_estimate_ud_spatstat(iteration = iteration, 
                            extract_coord = function(s) s$states,
                            cl = NULL, 
                            plot = FALSE)
# (optional) Examine selected UDs
lapply_qplot_ud(iteration, "spatstat", "h", "ud.tif")


#### End of code.
###########################
###########################