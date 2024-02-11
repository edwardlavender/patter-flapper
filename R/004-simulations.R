###########################
###########################
#### simulations.R

#### Aims
# 1) Compare simulations with algorithm outputs on the bathymetry grid
# * This is used to understand particle degeneracy

#### Prerequisites
# 1) Process spatial datasets


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
library(patter)
library(prettyGraphics)
library(tictoc)

#### Load data 
dv::src()
bathy     <- terra::rast(here_data("spatial", "bathy.tif"))
mpa       <- readRDS(here_data("spatial", "mpa.rds"))

#### Global params
ssf()


###########################
###########################
#### Simulate data 

#### (optional) Test
test <- FALSE
if (test) {
  n_path <- 2L
  cl     <- 2L
} else {
  n_path <- 20L
  cl     <- 20L
}

#### Define study area
# We zoom into the MPA for speed (e.g., for preparation of detection kernels)
boundary <- terra::ext(mpa)
bathy    <- terra::crop(bathy, boundary)
terra::plot(bathy)
terra::lines(mpa)

#### Define study period
n_step <- 60 * 24 * 2/2
period <- 
  seq(as.POSIXct("2016-01-01 09:00:00", tz = "UTC"),
      by = "2 mins", 
      length.out = n_step)

#### Define parameters
mobility <- 500
gamma    <- 750

#### Simulate paths (~3 mins, 20 paths)
# We will simulate N random walks (with different origins)
# We assume the individual was 'detected' at the start and end of each path
tic()
ssf()
paths <- 
  lapply(1:n_path, function(i) {
    tryCatch(
      {
        p <- sim_path_walk(.bathy = bathy, 
                           .origin = NULL, 
                           .lonlat = FALSE, 
                           .n_step = length(period), 
                           .n_path = 1L, 
                           .plot = FALSE)
        path_id <- timestamp <- NULL
        p[, path_id := i]
        p[, timestamp := period]
        # Validate mobility on the grid
        stopifnot(all(dist_along_path(cbind(p$cell_x, p$cell_y)) <= mobility, na.rm = TRUE))
        p
      },
      error = function(e) NULL)
  }) |> 
  plyr::compact()
toc()
qs::qsave(paths, here_data("sims", "paths.qs"))
# Plot an example path
head(paths[[1]])
n_path <- length(paths)
terra::plot(bathy)
add_sp_path(paths[[1]]$cell_x, paths[[1]]$cell_y, length = 0.01)

#### Define acoustic array(s)
arrays <- 
  lapply(paths, function(d) {
    data.table(receiver_id = c(1L, 2L), 
               receiver_easting = c(d$cell_x[1], d$cell_x[n_step]),
               receiver_northing = c(d$cell_y[1], d$cell_y[n_step]), 
               receiver_start = as.Date(min(period)), 
               receiver_end = as.Date(max(period)) + 1, 
               receiver_range = gamma)
  })

#### Simulate acoustic observations
acoustics <- 
  lapply(paths, function(d) {
    data.table(receiver_id = c(1L, 2L), 
               timestamp = d$timestamp[c(1L, n_step)])
  })
qs::qsave(acoustics, here_data("sims", "acoustics.qs"))

#### Simulate archival observations
archivals <- 
  lapply(paths, function(d) {
    ssf()
    d |> 
      select(timestamp, cell_z) |> 
      mutate(depth = runif(n(), cell_z - 5, cell_z + 5), 
             depth = if_else(depth < 0, 0, depth)) |> 
      as.data.table()
  })
qs::qsave(archivals, here_data("sims", "archivals.qs"))

#### Prepare algorithms (~6 mins, 20 paths)
tic()
dinputs <- 
  lapply(seq_len(n_path), function(i) {
  
  # Prepare data list
  dlist <- pat_setup_data(.acoustics = acoustics[[i]], 
                          .moorings = arrays[[i]], 
                          .archival = archivals[[i]], 
                          .bathy = bathy, 
                          .lonlat = FALSE)
  # Prepare observation timeline 
  obs <- pf_setup_obs(.dlist = dlist, 
                      .step = "2 mins",
                      .mobility = mobility, 
                      .receiver_range = gamma)
  obs[, depth_shallow_eps := 5]
  obs[, depth_deep_eps := 5]
  # Origin
  shallow <- dlist$spatial$bathy - obs$depth_shallow_eps[1]
  deep    <- dlist$spatial$bathy + obs$depth_deep_eps[1]
  origin  <- (obs$depth[1] >= shallow) & (obs$depth[1] <= deep)
  origin  <- terra::classify(origin, cbind(0, NA))
  if (FALSE) {
    terra::plot(origin)
    points(paths[[i]]$x[1], paths[[i]]$y[1])
  }
  dlist$spatial$origin <- origin
  # AC* likelihood components
  dlist$algorithm$detection_overlaps <- acs_setup_detection_overlaps(.dlist = dlist)
  dlist$algorithm$detection_kernels  <- acs_setup_detection_kernels(.dlist = dlist)
  
  # Outputs
  list(dlist = dlist, obs = obs)
  
})
toc()
# Check depths
sapply(paths, \(d) range(d$cell_z))


###########################
###########################
#### Run algorithms

#### Timings
# 20 paths, 20 cl, 20 mins

#### Define simulation run & inputs
tic()
success <- cl_lapply(seq_len(n_path), .cl = cl, .fun = function(i) {
# i <- 1L
message(paste0(rep("-", 25), collapse = ""))
print(i)
path     <- paths[[i]]
moorings <- arrays[[i]]
obs      <- dinputs[[i]]$obs
dlist    <- dinputs[[i]]$dlist

#### Forward run
# Define output connections
sink_pff <- here_data("sims", i, "forward")
unlink(sink_pff, recursive = TRUE)
dir.create(sink_pff, recursive = TRUE)
# Run simulation 
print("Running forward simulation...")
tic()
ssf()
out_pff <- 
  pf_forward(.obs = obs, 
             .dlist = dlist, 
             .likelihood = list(pf_lik_dc = pf_lik_dc, 
                                acs_filter_container = acs_filter_container), 
             .n = 1e3L, 
             .trial = pf_opt_trial(.trial_sampler = TRUE, 
                                   .trial_sampler_crit = 10L,
                                   .trial_revert = 0L), 
             .control = pf_opt_control(.sampler_batch_size = 1e3L),
             .record = pf_opt_record(.save = TRUE), 
             .verbose = file.path(sink_pff, "log.txt"))
toc()
# Write outputs to file (~0.2 s)
tic()
qs::qsave(out_pff, file.path(sink_pff, "out_pff.qs"))
pf_files_size(sink_pff)
toc()
# Validate convergence
data.frame(n_step_algorithm = length(out_pff$history), n_step = n_step)

#### Investigation
# * Can we achieve convergence by boosting n? 

#### Animations
# Examine region
terra::plot(bathy)
add_sp_path(path$cell_x, path$cell_y, length = 0.01)
# Define region
pext <- 
  out_pff$history |> 
  rbindlist(fill = TRUE) |> 
  summarise(xmin = min(c(x_now, path$x)), xmax = max(c(x_now, path$x)), 
            ymin = min(c(y_now, path$y)), ymax = max(c(y_now, path$y))) |> 
  as.data.table()
fig_dlist <- dlist
fig_dlist$spatial$bathy <- terra::crop(fig_dlist$spatial$bathy, unlist(pext))
# Make plots (~1.5 mins)
tic()
sink_fig <- here_fig("sims", i, "frames")
unlink(sink_fig, recursive = TRUE)
dir.create(sink_fig, recursive = TRUE)
tic()
print("Plotting frames...")
pf_plot_history(.dlist = fig_dlist,
                .forward = out_pff,
                .steps = NULL,
                .png = list(filename = sink_fig, height = 5, width = 5, res = 300),
                .add_forward = list(pch = 21, col = "red", bg = "red", cex = 0.25),
                .add_layer = function(t) {
                  # Plot true path
                  add_sp_path(path$x, path$y, length = 0.01, lwd = 0.75, 
                              col = scales::alpha("dimgrey", 0.75))
                  # Add moorings
                  points(moorings$receiver_easting, moorings$receiver_northing,
                         col = "blue", lwd = 2)
                  # Plot true position at time t
                  points(path$x[t], path$y[t], col = "black", lwd = 3)
                  terra::sbar(d = mobility, xy = "bottomleft")
                },
                .cl = 10L)
toc()
# Make animation (~1 min)
print("Animating frames...")
tic()
sink_mp4 <- here_fig("sims", i, "mp4")
unlink(sink_mp4, recursive = TRUE)
dir.create(sink_mp4, recursive = TRUE)
input  <- unlist(pf_files(sink_fig))
output <- file.path(sink_mp4, "ani.mp4")
av::av_encode_video(input, output)
toc()

#### Outputs
out_pff$convergence

})
toc()

#### Examine success
# Fails:
# * 2, 5, 6, 7, 13, 14
success <- unlist(success)
names(success) <- seq_len(n_path)
table(success)
success


#### End of code.
###########################
###########################