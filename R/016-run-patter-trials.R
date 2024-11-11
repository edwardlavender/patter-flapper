###########################
###########################
#### run-patter-trials.R

#### Aims
# 1) Trial convergence solutions for ACDC algorithm runs
# > Convergence is hard to achieve with moderate numbers of particles with the ACDC algorithm
# > Here we explore various options to achieve convergence at sensible computation times:

# * Modify tuning parameters
# - Boost particles 
# - Reduce resampling 
# * Modify movement model
# - Use simple random walk model 
# * Modify likelihood
# - Aggregate bathymetry & boost uncertainty in Normal depth model 
# - Remove (inflate) deep truncation parameter in Normal depth model
# - Implementation of depth data when resting only
# - Use deep truncation parameter and uniform depth model

#### Prerequisites
# 1) Previous scripts


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
# try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(truncdist)
dv::src()

#### Load data 
bb                <- qreadext(here_data("spatial", "bb.qs"))
iteration         <- qs::qread(here_data("input", "iteration", "patter.qs"))
moorings          <- qs::qread(here_data("input", "mefs", "moorings.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))
archival_by_unit  <- qs::qread(here_data("input", "archival_by_unit.qs"))
behaviour_by_unit <- qs::qread(here_data("input", "behaviour_by_unit.qs"))
mpa               <- qreadvect(here_data("spatial", "mpa.qs"))

#### Julia connect
julia_connect()
set_seed()
julia_command('Rasters.checkmem!(false);')
# TO DO
# Incorporate ModelObsCaptureContainer in workflow below
set_ModelObsCaptureContainer()
set_model_move_components()



###########################
###########################
#### Trials

#### Define iterations
# Define ACDC iterations 
iteration <- iteration[sensitivity == "best" & dataset == "acdc", ]
# (optional) Select individuals with 'good' time series which probably remained near MPA
# iteration <- iteration[unit_id %in% c(47, 49, 52, 74, 119), ] 
# (optional) Focus on individuals with convergence challenges 
# failures <- c(1, 20, 29, 30, 46, 57, 58, 72, 74, 109, 110, 115, 118, 121, 124)
# iteration <- iteration[unit_id %in% failures, ]
iteration

#### Define map 
map  <- terra::rast(here_data("spatial", "bathy-5m.tif"))
if (TRUE) {
  
  # Visualise map, incl. spikes
  # terra::plot(map)
  # terra::plot(map > 100)
  # terra::sbar(500)
  
  # Aggregate bathy onto a ~500 x 500 grid
  rnow <- terra::res(map)
  rnew <- 
    terra::rast(terra::ext(map), nrow = 500, ncol = 500, crs = terra::crs(map)) |> 
    terra::res()
  rnew <- c(rnew[2], rnew[1])
  map_agg <- terra::aggregate(map, fact = rnew / rnow, fun = "mean", na.rm = TRUE)
  map_agg
  
  # Visualise aggregation
  # > There are very few spikes around the coastline (good)
  terra::plot(map_agg > 100)
  
  # Change range within aggregated cells
  map_agg_max <-  terra::aggregate(map, fact = rnew / rnow, fun = "max", na.rm = TRUE)
  map_agg_min <- terra::aggregate(map, fact = rnew / rnow, fun = "min", na.rm = TRUE)
  map_agg_rng <- map_agg_max - map_agg_min
  terra::hist(map_agg_rng)
  terra::global(map_agg_rng, "sd", na.rm = TRUE)
  terra::global(map_agg_rng, quantile, probs = seq(0, 1, by = 0.01), na.rm = TRUE)
  
  # Check SD within aggregated cells 
  map_agg_sd <- terra::aggregate(map, fact = rnew / rnow, fun = "sd", na.rm = TRUE)
  terra::global(map_agg_sd, "mean", na.rm = TRUE)
  terra::global(map_agg_sd, quantile, prob = seq(0.95, 1, by = 0.001), na.rm = TRUE)
  
  # Visualise SD
  # * This is high along coastlines (steep relief)
  terra::plot(map_agg_sd)
  
  # Check maximum depth below depth of aggregated cell (+ 20 m)
  # * The differences are greatest near coastline
  map_agg_max   <- terra::aggregate(map, fact = rnew / rnow, fun = "max", na.rm = TRUE)
  map_agg_delta <- map_agg_max - map_agg
  terra::plot(map_agg_delta)
  terra::global(map_agg_delta, quantile, prob = seq(0.9, 1, by = 0.001), na.rm = TRUE) # 43.20033 m + 20 m
  
  # Use aggregated map
  map <- map_agg
  
}
names(map) <- "map_value"
map

#### Visualise maps
# Visualise study area
terra::plot(map)
points(moorings$receiver_x, moorings$receiver_y)
terra::sbar(20000) # 20 km
# Consider restricted region in which to sample initial locations
# bbinit <- matrix(c(681197.1, 6180218,
#                    681197.1, 6280651,
#                   723437.9, 6280651,
#                   723437.9, 6180218), 
#                  nrow = 4, byrow = TRUE) 
# lines(bbinit)
# origin <- terra::crop(map, terra::ext(bbinit))
# terra::plot(origin)

#### Set maps 
# (optional) Set restricted origin
# set_map(origin, .as_Raster = TRUE, .as_GeoArray = FALSE)
# Set wider maps
# set_map(map, .as_Raster = FALSE, .as_GeoArray = TRUE)
set_map(map)
# set_vmap(.map = map, .mobility = 1095, .plot = TRUE)

#### Iterate over example individuals and check convergence
tic()
for (i in 1:nrow(iteration)) {

  sim <- iteration[i, ]
  
  print(paste0(rep("-", 50), collapse = ""))
  print(paste(i, ":", sim$unit_id))
  
  #### Define individual datasets

  detections  <- acoustics_by_unit[[sim$unit_id]]
  archival    <- archival_by_unit[[sim$unit_id]]
  behaviour   <- behaviour_by_unit[[sim$unit_id]]
  # xinit     <- NULL
  xinit       <- qs::qread(here_data("input", "xinit", rho, sim$individual_id, sim$month_id, "xinit-fwd.qs"))
  # (optional) Trial implementation of depth data only when resting
  # archival <- archival[which(behaviour == 1L), ]
  
  #### Define timeline
  timeline <- patter_timeline(sim$month_id)
  stopifnot(length(timeline) == length(behaviour))
  
  #### Define movement model 
  state      <- state_flapper
  model_move <- move_flapper()
  julia_assign("behaviour", behaviour)
  update_model_move_components()
  
  # Visualise movement model realisations 
  print(rho)
  if (FALSE) {
    paths <- sim_path_walk(.map = map, 
                           .timeline = timeline, 
                           .state = state, 
                           .model_move = model_move,
                           .n_path = 4L, .one_page = TRUE)
    # Check correlation coefficient
    paths |> 
      group_by(path_id) |> 
      summarise(rho = circular::cor.circular(angle, dplyr::lead(angle))) |> 
      as.data.table() |> 
      suppressWarnings()
  }

  #### Define acoustic observations 
  moorings[, receiver_alpha := sim$receiver_alpha]
  moorings[, receiver_beta := sim$receiver_beta]
  moorings[, receiver_gamma := 3000]
  acoustics <- assemble_acoustics(.timeline = timeline,
                                  .detections = detections, 
                                  .moorings = moorings)
  
  #### Define acoustic containers
  containers <- assemble_acoustics_containers(.timeline  = timeline, 
                                              .acoustics = acoustics,
                                              .mobility = sim$mobility, 
                                              .map = map)
  
  #### Define archival observations
  archival <- assemble_archival(.timeline = timeline, 
                                .archival = 
                                  archival |> 
                                  rename(obs = "depth") |> 
                                  mutate(sensor_id = 1L, 
                                         depth_sigma = sim$depth_sigma, 
                                         depth_deep_eps = sim$depth_deep_eps) |>
                                  select("timestamp", "sensor_id", "obs", 
                                         "depth_sigma", "depth_deep_eps") |> 
                                  as.data.table())
  # Update parameters, e.g.:
  # * (20, 20): default
  # * (35, 20): inflated sigma (20 + 15 for 95 % quantile due to aggregation)
  # * (35, 500): inflated depth_deep_eps threshold due to aggregation 
  # * (100, 100): effectively 'uniform' probabilities above seabed depth reflecting unknown pelagic behaviour
  archival[, depth_sigma := 100]
  archival[, depth_deep_eps := 350]
  if (FALSE) {
    # seabed <- 0
    # seabed <- 100
    seabed <- 350
    curve(dtrunc(x, spec = "norm", a = 0, b = 350, mean = seabed, sd = 100), from = 0, to = 350)
  }

  #### Collect yobs
  # Define yobs for ACDC
  alg <- "acdc"
  yobs_fwd <- list(ModelObsAcousticLogisTrunc = acoustics,
                   ModelObsAcousticContainer = containers$forward, 
                   ModelObsDepthNormalTrunc = archival)
  yobs_bwd <- list(ModelObsAcousticLogisTrunc = acoustics,
                   ModelObsAcousticContainer = containers$backward, 
                   ModelObsDepthNormalTrunc = archival)
  # (optional) Define t_resample 
  t_resample_fwd <- sort(unique(c(which(timeline %in% containers$forward$timestamp), which(timeline %in% acoustics$timestamp[acoustics$obs == 1L]))))
  t_resample_bwd <- sort(unique(c(which(timeline %in% containers$backward$timestamp), which(timeline %in% acoustics$timestamp[acoustics$obs == 1L]))))
  if (length(t_resample_fwd) == 0L) {
    t_resample_fwd <- NULL
  }
  if (length(t_resample_bwd) == 0L) {
    t_resample_bwd <- NULL
  }
  
  # (optional) Set yobs for AC
  # alg <- "ac"
  if (alg == "ac") {
    yobs_fwd$ModelObsDepthNormalTrunc <- NULL
    yobs_bwd$ModelObsDepthNormalTrunc <- NULL
  }
  # (optional) Set yobs for DC
  alg <- "dc"
  if (alg == "dc") {
    t_resample_fwd <- t_resample_bwd <- NULL
    yobs_fwd$ModelObsAcousticLogisTrunc <- NULL
    yobs_fwd$ModelObsAcousticContainer  <- NULL
    yobs_bwd$ModelObsAcousticLogisTrunc <- NULL
    yobs_bwd$ModelObsAcousticContainer  <- NULL
  }
  
  #### Collect args 
  # n_resample     <- as.numeric(500)
  t_resample_fwd <- NULL
  n_resample     <- 1e5
  args <- list(.timeline   = timeline,
               .state      = state,
               .xinit      = xinit,
               .yobs       = yobs_fwd,
               .model_move = model_move,
               .n_particle = 5000L, 
               .n_move     = 10000L,
               .n_resample = n_resample,
               .t_resample = t_resample_fwd,
               .n_iter     = 1L,
               .direction  = "forward",
               .verbose    = TRUE)
  
  #### Run filter and check convergence
  set_seed()
  fwd <- tryCatch(do.call(pf_filter, args, quote = TRUE), 
                  error = function(e) e)
  if (inherits(fwd, "error")) {
    message(fwd$message)
  }
  
  #### (optional) Run downstream analyses
  # For an example individual run downstream analyses
  # Compare maps for AC, DC and ACDC
  # Visualise the extent to which the incorporation of depth parameters improve AC maps
  if (TRUE) {
    # Run backward filter 
    args$.xinit        <- qs::qread(here_data("input", "xinit", rho, sim$individual_id, sim$month_id, "xinit-fwd.qs"))
    args$.yobs         <- yobs_bwd
    # args$.t_resample <- t_resample_bwd
    args$.direction    <- "backward"
    set_seed()
    bwd <- do.call(pf_filter, args, quote = TRUE)
  }
  
  if (FALSE) {
    # Run smoother
    set_seed()
    smo <- pf_smoother_two_filter(.n_particle = 500L)
    
    # Visualise map for algorithm in fig/trials/
    png(here_fig("trials", glue("{sim$index}-{alg}.png")), 
        height = 10, width = 10, units = "in", res = 600)
    ud <- map_pou(.map = map, .coord = smo$states)$ud
    points(moorings$receiver_x, moorings$receiver_y)
    # map_hr_core(.map = ud, .add = TRUE)
    dev.off()
  }
  
  
}
toc()
beepr::beep(10)


#### End of code. 
###########################
###########################