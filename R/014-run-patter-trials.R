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
dv::src()

#### Load data 
bb                <- qreadext(here_data("spatial", "bb.qs"))
iteration         <- qs::qread(here_data("input", "iteration", "patter.qs"))
moorings          <- qs::qread(here_data("input", "mefs", "moorings.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))
archival_by_unit  <- qs::qread(here_data("input", "archival_by_unit.qs"))
behaviour_by_unit <- qs::qread(here_data("input", "behaviour_by_unit.qs"))

#### Julia connect
julia_connect()
set_seed()
julia_command(ModelMoveFlapper)


###########################
###########################
#### Trials

#### Define ACDC iterations
iteration <- iteration[sensitivity == "best" & dataset == "acdc", ]

#### Define map 
map  <- terra::rast(here_data("spatial", "bathy-5m.tif"))
if (TRUE) {
  rnow <- terra::res(map)
  rnew <- 
    terra::rast(bb, nrow = 500, ncol = 500, crs = terra::crs(map)) |> 
    terra::res()
  map_agg <- terra::aggregate(map, fact = rnew / rnow, fun = "mean", na.rm = TRUE)
  terra::plot(map_agg)
  map_agg_sd <- terra::aggregate(map, fact = rnew / rnow, fun = "sd", na.rm = TRUE)
  terra::global(map_agg_sd, "mean")
  map <- map_agg
}
names(map) <- "map_value"
set_map(map)
map

#### Iterate over example individuals and check convergence
for (i in 1:5) {
  
  print(paste0(rep("-", 50), collapse = ""))
  print(i)
  
  #### Define individual datasets
  sim         <- iteration[i, ]
  detections  <- acoustics_by_unit[[sim$unit_id]]
  archival    <- archival_by_unit[[sim$unit_id]]
  behaviour   <- behaviour_by_unit[[sim$unit_id]]
  # (optional) Trial implementation of depth data only when resting
  # archival <- archival[which(behaviour == 1L), ]
  
  #### Define timeline
  timeline <- patter_timeline(sim$month_id)
  
  #### Define movement model 
  state       <- "StateXY"
  model_move  <- patter_ModelMove(sim)
  julia_assign("behaviour", behaviour)
  JuliaCall::julia_command(simulate_step.ModelMoveFlapper)
  JuliaCall::julia_command(logpdf_step.ModelMoveFlapper)
  
  #### Define acoustic observations 
  moorings[, receiver_alpha := sim$receiver_alpha]
  moorings[, receiver_beta := sim$receiver_beta]
  moorings[, receiver_gamma := sim$receiver_gamma]
  acoustics <- assemble_acoustics(.timeline = timeline,
                                  .detections = detections, 
                                  .moorings = moorings)
  
  #### Define acoustic containers
  containers <- assemble_acoustics_containers(.timeline  = timeline, 
                                              .acoustics = acoustics,
                                              .mobility = sim$mobility, 
                                              .threshold = 139199)
  
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
  archival[, depth_sigma := 20]
  archival[, depth_deep_eps := 20]
  curve(dnorm(x, 0, 40), from = 0, to = 350)
  curve(dnorm(x, 200, 40), from = 0, to = 350)
  
  #### Collect yobs
  yobs <- list(ModelObsAcousticLogisTrunc = acoustics,
               ModelObsAcousticContainer = containers$forward, 
               ModelObsDepthNormalTrunc = archival)
  
  #### Collect args 
  args <- list(.timeline   = timeline,
               .state      = state,
               .xinit      = NULL,
               .yobs       = yobs,
               .model_move = model_move,
               .n_particle = 1e5L, 
               .n_move     = 10000L,
               .n_iter     = 1L,
               .direction  = "forward",
               .verbose    = TRUE)
  
  #### Run filter and check convergence
  pout <- do.call(pf_filter, args, quote = TRUE) 
  
}

#### Results (aggregated bathymetry)
# * Normal depth model (20, 20): 
# * Normal depth model (40, 20):
# * Normal depth model (20, 200):
# * Normal depth model (20, 200):
# * Uniform depth model (-100, 0):
# * Uniform depth model (-100, 0):


#### End of code. 
###########################
###########################