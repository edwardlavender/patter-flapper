###########################
###########################
#### simulate-algorithms.R

#### Aims
# 1) Simulation exploration & validation

#### Prerequisites
# 1) NA

#### Work flow
#
# * Simulate three movement paths 
# * Estimate the 'true' UD
# * Simulate three corresponding acoustic & archival observational datasets
#
# * Run the COA algorithm
# - Use three delta t values (based on trial & error selection)
#
# * Run the RSP algorithm
# - Use three er.add values (based on trial & error selection)
#
# * Run the patter algorithms
# - For each dataset/algorithm (ACPF, DCPF, ACDCPF)
# - For each set of the number of particles (10,000, 20,000, ... etc.)
# - For each set of parameter values
# - For 1, ..., 30 iterations
# --> Run filter & record convergence (for convergence analysis)
# --> If also maximum number of particles, run smoother (for residency analysis)
# --> If also iterations 1:3, estimate UDs (for maps)
# --> For UD estimation, use coarse npixel for speed


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
map      <- terra::rast(here_data("spatial", "bathy.tif"))
win      <- qs::qread(dv::here_data("spatial", "win.qs"))
skateids <- qs::qread(here_data("input", "mefs", "skateids.qs"))
# moorings <- qs::qread(here_data("input", "mefs", "moorings.qs"))
# Best-guess ModelObs parameters
# * Define in formulate-models.R
moorings <- qs::qread(here_data("input", "mefs", "moorings.qs")) 
ModelObsAcousticLogisTruncPars <- 
  moorings |> 
  select(sensor_id = "receiver_id", "receiver_x", "receiver_y") |> 
  mutate(receiver_alpha = 4, receiver_beta = -0.01, receiver_gamma = 1750) |> 
  as.data.table()
ModelObsDepthNormalTruncPars <- data.table(sensor_id = 1, 
                                           sigma = 20, 
                                           depth_deep_eps = 20)
ModelObsPars <- list(ModelObsAcousticLogisTrunc = ModelObsAcousticLogisTruncPars, 
                     ModelObsDepthNormalTrunc = ModelObsDepthNormalTruncPars)

#### Julia set up
julia_connect()
set_seed()
set_map(map)
julia_command(ModelObsAcousticContainer)
julia_command(ModelObsAcousticContainer.logpdf_obs)


###########################
###########################
#### Set up simulations 

#### Set up datasets
# Define tagging locations in UTM 29N
xinits <- 
  cbind(skateids$long_tag_capture, skateids$lat_tag_capture) |>
  terra::vect(crs = "EPSG:4326") |>
  terra::project(terra::crs(map)) |> 
  terra::crds() |>
  as.data.table() |>
  mutate(map_value = terra::extract(map, cbind(x, y))[, 1]) |> 
  select("map_value", "x", "y") |> 
  as.data.table()

#### COA algorithm
# Define parameters
pcoa <- data.table(parameter_id = 1:3L, 
                   delta_t = c("6 hours", "12 hours", "24 hours"))
# Define iteration dataset
iteration_coa <- 
  CJ(path = 1:3L, parameter_id = 1:3, dataset = "ac", iter = 1L) |>
  mutate(dataset = "ac", 
         delta_t = pcoa$delta_t[match(parameter_id, pcoa$parameter_id)], 
         folder = file.path("data", "output", "simulation", "coa"),
         folder_coord = file.path(folder, dataset, parameter_id, iter, "coord"), 
         folder_ud = file.path(folder, dataset, parameter_id, iter, "ud")
         ) |> 
  as.data.table()
# Build folders
# dirs.create(iteration_coa$folder_coord)
# dirs.create(iteration_coa$folder_ud)
# dirs.create(file.path(iteration_coa$folder_ud, "spatstat", "h"))

#### RSP algorithm

#### Patter algorithms 

#### Estimation parameters

spatstat.geom::spatstat.options("npixel" = 100)
 



###########################
###########################
#### Simulate paths and observations

#### Define simulation timeline 
timeline <- seq(as.POSIXct("2024-04-01 00:00:00", tz = "UTC"), 
                as.POSIXct("2024-04-30 23:58:00", tz = "UTC"), 
                by = "2 mins")

#### Simulate initial location
# We could sample from the tagging locations
# But we instead sample from the receiver array to ensure  
# ... we generate some detections for the COA/RSP/AC*PF algorithms
# xinit   <- xinits[sample.int(nrow(xinits), 1), ]
xinit_bb  <- terra::ext(min(moorings$receiver_x), max(moorings$receiver_x), 
                        min(moorings$receiver_y), max(moorings$receiver_y))
xinit_map <- terra::crop(map, xinit_bb)
xinit     <- terra::spatSample(xinit_map, size = 1L, xy = TRUE, na.rm = TRUE)
xinit     <- data.table(map_value = xinit$map_value, x = xinit$x, y = xinit$y)


#### (Temporary) resting investigation
archival <- qs::qread(here_data("input", "archival_by_unit.qs")) |> plyr::compact()
archival <- rbindlist(archival)
archival[, fct := paste(individual_id, mmyy)]
archival[, state := get_mvt_resting(copy(archival), fct = "fct")]


arc <- archival[1:1000, ]

duration <- rle(arc$state)
duration <- duration$lengths[duration$values == 0L]

fitdistrplus::descdist(duration)

dbn <- "gamma"
dbn <- "exponential"
pars <- MASS::fitdistr(duration, "gamma")
plotdist(duration, "gamma", as.list(pars$estimate))

# hist(duration, breaks = 5, probability = TRUE)
# curve(dcauchy(x, pars$estimate[1], pars$estimate[2]), col = "red", lwd = 2, add = TRUE)

#### Define movement model
behaviour <- sample(c(1L, 2L), length(timeline), replace = TRUE)
julia_assign("behaviour", behaviour)
julia_command(ModelMoveFlapper)
julia_command(simulate_step.ModelMoveFlapper)
julia_command(logpdf_step.ModelMoveFlapper)



coord_path <- sim_path_walk(.map = map, 
                            .timeline = timeline, 
                            .state = "StateXY", 
                            .xinit = xinit, 
                            .model_move = move_flapper(), 
                            .n_path = 1L, 
                            .plot = TRUE)




points(moorings$receiver_x, moorings$receiver_y)
# Map path UD (~58 s with 100 pixels)
ud_path <- map_dens(.map = map, 
                    .owin = win,
                    sigma = bw.h,
                    .coord = coord_path[, .(x, y)]
)
# Simulate observations
yobs <- sim_observations(.timeline = timeline, 
                         .model_obs = c("ModelObsAcousticLogisTrunc", "ModelObsDepthNormalTrunc"), 
                         .model_obs_pars = ModelObsPars)

###########################
###########################
#### Run algorithms

#### COA algorithm
# Estimate coordinates 
table(yobs$ModelObsAcousticLogisTrunc[[1]]$obs)
coord_coa <- coa(.map = map, 
                 .acoustics = yobs$ModelObsAcousticLogisTrunc[[1]][obs > 0L, ],
                 .delta_t = "24 hours") 
stopifnot(nrow(coord_coa) > 0)
terra::plot(map)
points(coord_coa$x, coord_coa$y)
# Estimate UD
ud_path <- map_dens(.map = map, 
                    .owin = win,
                    sigma = bw.h,
                    .coord = coord_coa[, .(x, y)]
                    )
