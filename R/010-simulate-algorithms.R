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
# - For 1, ..., N iterations
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
library(spatstat)
library(ggplot2)

#### Load data 
map      <- terra::rast(here_data("spatial", "bathy.tif"))
ud_grid  <- terra::rast(here_data("spatial", "ud-grid.tif"))
ud_null  <- terra::rast(here_data("spatial", "ud-null.tif"))
win      <- qs::qread(dv::here_data("spatial", "win.qs"))
skateids <- qs::qread(here_data("input", "mefs", "skateids.qs"))
moorings <- qs::qread(here_data("input", "mefs", "moorings.qs")) 
pars     <- qs::qread(here_data("input", "pars.qs"))


###########################
###########################
#### Set up algorithms 

#### Julia set up
julia_connect()
set_seed()
set_map(map)
julia_command(ModelMoveFlapper)
julia_command(ModelObsAcousticContainer)
julia_command(ModelObsAcousticContainer.logpdf_obs)

#### Define simulation settings
# We simulate three paths
n_path   <- 3L
# Simulation timeline
timeline <- seq(as.POSIXct("2024-04-01 00:00:00", tz = "UTC"), 
                as.POSIXct("2024-04-30 23:58:00", tz = "UTC"), 
                by = "2 mins")

#### Define ModelMove structure to simulate paths
# We simulate the path using the best-guess parameters
model_move <- patter_ModelMove(pars$pmovement[1, ])

#### Define ModelObs structures to simulate observations 
# We simulate parameters using best-guess parameters
ModelObsAcousticLogisTruncPars <- 
  moorings |> 
  select(sensor_id = "receiver_id", "receiver_x", "receiver_y") |> 
  mutate(receiver_alpha = pars$pdetection$receiver_alpha[1], 
         receiver_beta = pars$pdetection$receiver_beta[1], 
         receiver_gamma = pars$pdetection$receiver_gamma[1]) |> 
  as.data.table()
ModelObsDepthNormalTruncPars <- data.table(sensor_id = 1L, 
                                           sigma = pars$pdepth$depth_sigma, 
                                           depth_deep_eps = pars$pdepth$depth_deep_eps)
ModelObsPars <- list(ModelObsAcousticLogisTrunc = ModelObsAcousticLogisTruncPars, 
                     ModelObsDepthNormalTrunc = ModelObsDepthNormalTruncPars)

#### Define ancillary datasets for algorithm runs
# COA & patter run with simulated acoustic data
# RSP requires moorings dataset to prepare the actel dataset
# * Receiver deployment dates must be dates that align with simulation timeline
moorings <- 
  moorings |> 
  mutate(receiver_start = as.Date(min(timeline) - lubridate::days(1)), 
         receiver_end = as.Date(as.Date(max(timeline) + lubridate::days(1))), 
         receiver_gamma = pars$pdetection$receiver_gamma[1]) |> 
  as.data.table()
head(moorings)


###########################
###########################
#### Simulate paths and observations

if (FALSE) {
  
  iteration_path <- data.table(s = c(2, 3, 4), id = c(1, 2, 3))
  lapply(split(iteration_path, 1:3), function(sim) {
    
    # This code is run three times, generating three paths & corresponding observational datasets
    # s <- 2; id <- 1 # path 1
    # s <- 3; id <- 2 # path 2
    # s <- 4; id <- 3 # path 3
    s <- sim$s; id <- sim$id
    set_seed(s)
    
    #### Simulate initial location
    # We sample an initial location from the receiver array 
    # This helps to ensure we generate some detections for the COA/RSP/AC*PF algorithms
    xinit_bb  <- terra::ext(min(moorings$receiver_x), max(moorings$receiver_x), 
                            min(moorings$receiver_y), max(moorings$receiver_y))
    xinit_map <- terra::crop(map, xinit_bb)
    xinit     <- terra::spatSample(xinit_map, size = 1L, xy = TRUE, na.rm = TRUE)
    xinit     <- data.table(map_value = xinit$map_value, x = xinit$x, y = xinit$y)
    
    #### Simulate behavioural states 
    # behaviour <- sample(c(1L, 2L), length(timeline), replace = TRUE)
    behaviour   <- simulate_behaviour(timeline)
    julia_assign("behaviour", behaviour)
    
    #### Define movement model
    julia_command(simulate_step.ModelMoveFlapper)
    julia_command(logpdf_step.ModelMoveFlapper)
    
    #### Simulate movement path 
    coord_path <- sim_path_walk(.map = map, 
                                .timeline = timeline, 
                                .state = "StateXY", 
                                .xinit = xinit, 
                                .model_move = model_move, 
                                .n_path = 1L, 
                                .plot = TRUE)
    points(moorings$receiver_x, moorings$receiver_y)
    
    #### Simulate observations
    yobs <- sim_observations(.timeline = timeline, 
                             .model_obs = c("ModelObsAcousticLogisTrunc", "ModelObsDepthNormalTrunc"), 
                             .model_obs_pars = ModelObsPars)
    # Check we have simulated detections
    stopifnot(length(which(yobs$ModelObsAcousticLogisTrunc[[1]]$obs == 1L)) > 100)
    table(yobs$ModelObsAcousticLogisTrunc[[1]]$obs)
    
    #### Map path UD 
    # * 13 s 500 pixels 
    stopifnot(spatstat.geom::spatstat.options("npixel") == 500)
    ud_path <- map_dens(.map = ud_grid, 
                        .owin = win,
                        .coord = coord_path[, .(x, y)],
                        .discretise = TRUE,
                        sigma = bw.h, 
                        .fterra = TRUE)
    
    #### Save datasets
    qs::qsave(behaviour, here_data("input", "simulation", id, "behaviour.qs"))
    qs::qsave(coord_path, here_data("input", "simulation", id, "coord.qs"))
    qs::qsave(yobs, here_data("input", "simulation", id, "yobs.qs"))
    terra::writeRaster(ud_path$ud, here_data("input", "simulation", id, "ud.tif"), 
                       overwrite = TRUE)
    beepr::beep(10)
    
  })
  
}

#### Load datasets
# Behaviour
behaviour_by_unit <- lapply(seq_len(n_path), function(path) {
  qs::qread(here_data("input", "simulation", path, "behaviour.qs"))
})
# Detections
detections_by_unit <- lapply(seq_len(n_path), function(path) {
  acoustics  <- qs::qread(here_data("input", "simulation", path, "yobs.qs"))
  acoustics$ModelObsAcousticLogisTrunc[[1]] |> 
    filter(obs == 1L) |>
    select(timestamp, receiver_id = sensor_id, receiver_x, receiver_y) |>
    as.data.table()
})
# Depths
archival_by_unit <- lapply(seq_len(n_path), function(path) {
  yobs <- qs::qread(here_data("input", "simulation", path, "yobs.qs"))
  yobs$ModelObsDepthNormalTrunc[[1]]
})


###########################
###########################
#### Run algorithms

if (FALSE) {
  
  #### Build iteration dataset
  # Define parameters
  pcoa <- data.table(parameter_id = 1:10L, 
                     delta_t = c("1 hour", "2 hours", "6 hours", "12 hours", "1 day", 
                                 "2 days", "3 days", "4 days", "5 days", "6 days"))
  # Define dataset
  iteration_coa <- 
    CJ(unit_id = seq_len(n_path), parameter_id = pcoa$parameter_id, dataset = "ac", iter = 1L) |>
    mutate(dataset = "ac", 
           delta_t = pcoa$delta_t[match(parameter_id, pcoa$parameter_id)], 
           folder = file.path("data", "output", "simulation", unit_id, "coa"),
           folder_coord = file.path(folder, dataset, parameter_id, iter, "coord"), 
           folder_ud = file.path(folder, dataset, parameter_id, iter, "ud")
    ) |> 
    arrange(unit_id, parameter_id) |>
    as.data.table()
  # Build folders
  dirs.create(iteration_coa$folder_coord)
  dirs.create(iteration_coa$folder_ud)
  dirs.create(file.path(iteration_coa$folder_ud, "spatstat", "h"))
  
  #### Define datasets
  datasets <- list(detections_by_unit = detections_by_unit, moorings = NULL)
  
  #### Run algorithm
  if (FALSE) {
    
    #### Estimate coordinates (~2 s)
    iteration <- copy(iteration_coa)
    iteration[, file_coord := file.path(folder_coord, "coord.qs")]
    lapply_estimate_coord_coa(iteration = iteration, datasets = datasets)
    # (optional) Examine selected coords
    lapply_qplot_coord(iteration, "coord.qs")
    
    #### Estimate UDs
    # Implementation (44 s, 500 pixels, sigma = bw.h, cl = 10L)
    nrow(iteration)
    lapply_estimate_ud_spatstat(iteration = iteration, 
                                extract_coord = NULL,
                                cl = 10L, 
                                plot = FALSE)
    # (optional) Examine selected UDs
    lapply_qplot_ud(iteration, "spatstat", "h", "ud.tif")
    
  }
  
  #### Visualise ME
  # Compute ME 
  me <- pbapply::pbsapply(split(iteration, seq_row(iteration)), function(sim) {
    # sim <- iteration[1, ]
    skill_me(.obs = terra::rast(here_data("input", "simulation", sim$unit_id, "ud.tif")), 
             .mod = terra::rast(file.path(sim$folder_ud, "spatstat", "h", "ud.tif")))
  })
  me_null <- pbapply::pbsapply(split(iteration, seq_row(iteration)), function(sim) {
    skill_me(.obs = terra::rast(here_data("input", "simulation", sim$unit_id, "ud.tif")), 
             .mod = ud_null)
  })
  # Visualise ME ~ delta_t
  it_me <- 
    iteration |>
    mutate(me = me, 
           delta_t = factor(delta_t, levels = pcoa$delta_t)) |> 
    as.data.table()
  it_me_avg <- 
    it_me |> 
    group_by(delta_t) |> 
    summarise(me = mean(.data$me)) |> 
    as.data.table()
  it_me |> 
    as.data.table() |> 
    ggplot(aes(delta_t, me, colour = factor(unit_id), group = unit_id)) + 
    geom_point() + 
    geom_line() + 
    geom_point(aes(er.ad, me_null)) + 
    geom_line(aes(er.ad, me_null)) + 
    geom_line(data = it_me_avg, aes(delta_t, me),
              colour = "black", group = 1) 
  
  # > Best guess:       3 days
  # > Restricted value: 2 day
  # > Flexible value:   4 days 
  
  ### Visualise maps
  # Define panel row (path) and column (parameter) labels
  cols <- c("NA", "3 days", "2 days", "4 days")
  cols <- factor(cols, levels = cols)
  rows <- c("Path", "COA[1]", "COA[2]", "COA[3]")
  rows <- factor(rows, levels = rows)
  panel <- data.table(row = rows, column = cols)
  # Define mapfiles for selected algorithm runs 
  mapfiles_alg <- 
    iteration |> 
    mutate(row = unit_id) |>
    filter(delta_t %in% panel$column) |> 
    mutate(column = panel$column[match(delta_t, panel$column)]) |>
    mutate(mapfile = file.path(folder_ud, "spatstat", "h", "ud.tif")) |>
    select(row, column, mapfile) |> 
    as.data.table()
  # Define mapfiles for simulated paths
  mapfiles_path <- 
    data.table(row = seq_len(n_path), 
               column = cols[1],
               mapfile = here_data("input", "simulation", seq_len(n_path), "ud.tif")
    )
  # Collect mapfiles (row, column, mapfile)
  mapfiles <- 
    rbind(mapfiles_alg, mapfiles_path) |> 
    arrange(row, column) |> 
    as.data.table()
  # Quick visual check of selected map 
  terra::plot(terra::rast(mapfiles$mapfile[5]))
  # Make maps (~5 s)
  ggplot_maps(mapdt = mapfiles, 
              png_args = list(filename = here_fig("simulation", "map-coa.png"), 
                              height = 5, width = 10, units = "in", res = 600))
  
}


###########################
###########################
#### RSP algorithms

if (FALSE) {
  
  #### Build iteration dataset
  # Define parameters
  prsp <- data.table(parameter_id = 1:22L, 
                     er.ad = c(250 * 0.01, 
                               250 * 0.05, # default 
                               250 * seq(0.1, 2, by = 0.1))
                     )
  # Define dataset
  iteration_rsp <- 
    CJ(unit_id = seq_len(n_path), parameter_id = prsp$parameter_id, dataset = "ac", iter = 1L) |>
    mutate(dataset = "ac", 
           er.ad = prsp$er.ad[match(parameter_id, prsp$parameter_id)], 
           folder = file.path("data", "output", "simulation", unit_id, "rsp"),
           folder_coord = file.path(folder, dataset, parameter_id, iter, "coord"), 
           folder_ud = file.path(folder, dataset, parameter_id, iter, "ud")
    ) |> 
    arrange(unit_id, parameter_id) |>
    as.data.table()
  # Build folders
  dirs.create(iteration_rsp$folder_coord)
  dirs.create(iteration_rsp$folder_ud)
  dirs.create(file.path(iteration_rsp$folder_ud, "spatstat", "h"))
  dirs.create(file.path(iteration_rsp$folder_ud, "dbbmm"))
  
  #### Define datasets
  datasets <- list(detections_by_unit = detections_by_unit, moorings = copy(moorings))
  
  #### Run algorithm
  if (FALSE) {
    
    #### Estimate coordinates (~2-4 mins)
    iteration <- copy(iteration_rsp)
    iteration[, file_coord := file.path(folder_coord, "coord.qs")]
    lapply_estimate_coord_rsp(iteration = iteration, datasets = datasets)
    # (optional) Examine selected coords
    lapply_qplot_coord(iteration, 
                       "coord.qs",
                       extract_coord = function(coord) {
                         cbind(coord$detections[[1]]$Longitude, 
                               coord$detections[[1]]$Latitude) |> 
                           terra::vect(crs = "EPSG:4326") |> 
                           terra::project("EPSG:32629") |> 
                           terra::crds() |> 
                           as.data.frame()
                       })
    
    #### Estimate UDs
    # Implementation (0.5-1.5 mins, 10 cl)
    lapply_estimate_ud_dbbmm(iteration = iteration, 
                             cl = 10L, 
                             plot = FALSE)
    # (optional) Examine selected UDs
    lapply_qplot_ud(iteration, "dbbmm", "ud.tif")
    
  }
  
  #### Visualise ME 
  me <- pbapply::pbsapply(split(iteration, seq_row(iteration)), function(sim) {
    # sim <- iteration[10, ]
    tryCatch(   
      skill_me(.obs = terra::rast(here_data("input", "simulation", sim$unit_id, "ud.tif")), 
               .mod = terra::rast(file.path(sim$folder_ud, "dbbmm", "ud.tif"))), 
      error = function(e) NA)
  })
  me_null <- pbapply::pbsapply(split(iteration, seq_row(iteration)), function(sim) {
      skill_me(.obs = terra::rast(here_data("input", "simulation", sim$unit_id, "ud.tif")), 
               .mod = ud_null)
  })
  # Visualise ME ~ er.ad
  it_me <- 
    iteration |>
    mutate(me = me, 
           me_null = me_null,
           er.ad = factor(er.ad, levels = prsp$er.ad)) |> 
    as.data.table()
  it_me_avg <- 
    it_me |> 
    group_by(er.ad) |> 
    summarise(me = mean(.data$me)) |> 
    as.data.table()
  it_me |> 
    as.data.table() |> 
    ggplot(aes(er.ad, me, colour = factor(unit_id), group = unit_id)) + 
    geom_point() + 
    geom_line() + 
    geom_point(aes(er.ad, me_null)) + 
    geom_line(aes(er.ad, me_null)) + 
    geom_line(data = it_me_avg, aes(er.ad, me),
              colour = "black", group = 1) 
  
  # > Best guess:       100
  # > Restricted value: 50
  # > Flexible value:   150 
  
  ### Visualise maps
  # Define panel row (path) and column (parameter) labels
  cols <- c("NA", "100", "25", "250")
  cols <- factor(cols, levels = cols)
  rows <- c("Path", "RSP[1]", "RSP[2]", "RSP[3]")
  rows <- factor(rows, levels = rows)
  panel <- data.table(row = rows, column = cols)
  # Define mapfiles for selected algorithm runs 
  mapfiles_alg <- 
    iteration |> 
    mutate(row = unit_id) |>
    filter(er.ad %in% panel$column) |> 
    mutate(column = panel$column[match(er.ad, panel$column)]) |>
    mutate(mapfile = file.path(folder_ud, "dbbmm", "ud.tif")) |>
    select(row, column, mapfile) |> 
    as.data.table()
  # Define mapfiles for simulated paths
  mapfiles_path <- 
    data.table(row = seq_len(n_path), 
               column = cols[1],
               mapfile = here_data("input", "simulation", seq_len(n_path), "ud.tif")
    )
  # Collect mapfiles (row, column, mapfile)
  mapfiles <- 
    rbind(mapfiles_alg, mapfiles_path) |> 
    arrange(row, column) |> 
    as.data.table()
  # Quick visual check of selected map 
  terra::plot(terra::rast(mapfiles$mapfile[5]))
  # Make maps
  ggplot_maps(mapdt = mapfiles, 
              png_args = list(filename = here_fig("simulation", "map-rsp.png"), 
                              height = 5, width = 10, units = "in", res = 600))
  
}


###########################
###########################
#### Patter algorithms

#### Define unitsets (unit_ids & algorithms)
unitsets <- 
  CJ(unit_id = seq_len(n_path), 
     dataset = c("ac", "dc", "acdc")) |>
  arrange(unit_id, factor(dataset, c("ac", "dc", "acdc"))) |>
  as.data.table()

#### Define parameters
# Observation model parameters 
pmovement  <- pars$pmovement
pdetection <- pars$pdetection
pdepth     <- pars$pdepth
parameters <- 
  rbind(
    cbind(sensitivity = "best", pmovement[1, ], pdetection[1, ], pdepth[1, ]),
    cbind(sensitivity = "movement", pmovement[2:3, ], pdetection[1, ], pdepth[1, ]),
    cbind(sensitivity = "ac", pmovement[1, ], pdetection[2:3, ], pdepth[1, ]),
    cbind(sensitivity = "dc", pmovement[1, ], pdetection[1, ], pdepth[2:3, ]) 
  )
parameters[, model_obs := seq_len(.N)]
# Algorithm settings (particles & no. of iters)
np <- c(10000, 25000, 50000, 100000, 250000)
ni <- 1:10L
# Collect algorithm settings and model parameters
parameters <- 
  CJ(model_obs = parameters$model_obs, np = np) |> 
  left_join(parameters, 
            by = "model_obs") |> 
  mutate(parameter_id = row_number()) |> 
  as.data.table()
# We run each setting ni times
parameters <- 
  lapply(ni, function(i) {
  p <- copy(parameters)
  p[, iter := i]
  p
}) |> rbindlist()

#### Define iteration dataset
iteration_patter <- lapply(split(unitsets, seq_len(nrow(unitsets))), function(d) {
  
  # Keep the relevant parameters, dependent upon the algorithm
  if (d$dataset == "ac") {
    p <- parameters[sensitivity %in% c("best", "movement", "ac"), ]
    p[, c("depth_sigma", "depth_deep_eps") := NA_real_]
  }
  if (d$dataset == "dc") {
    p <- parameters[sensitivity %in% c("best", "movement", "dc"), ]
    p[, c("receiver_alpha", "receiver_beta", "receiver_gamma") := NA_real_]
  }
  if (d$dataset == "acdc") {
    p <- copy(parameters)
  }
  cbind(d, p)
  
}) |> 
  rbindlist() |> 
  mutate(index = row_number(), 
         folder = file.path("data", "output", "simulation", unit_id, "patter"),
         folder_coord = file.path(folder, dataset, parameter_id, iter, "coord"), 
         folder_ud = file.path(folder, dataset, parameter_id, iter, "ud")) |>
  dplyr::select("index", 
         "unit_id",
         "dataset", 
         "parameter_id", "iter", "sensitivity", 
         "k1", "k2", "theta1", "theta2", "mobility", 
         "receiver_alpha", "receiver_beta",
         "receiver_gamma", "depth_sigma", "depth_deep_eps",
         "np",
         "folder_coord", "folder_ud") |>
  as.data.table()

#### Add supporting columns
# sim$month_id is required by estimate_coord_patter() to define the timeline
iteration_patter[, month_id := "042024"]
# Implement smoothing for simulations that use the max. number of particles
iteration_patter[, smooth := FALSE]
iteration_patter[np == max(np), smooth := TRUE]

#### Build patter folders
dirs.create(iteration_patter$folder_coord)
dirs.create(iteration_patter$folder_ud)
dirs.create(file.path(iteration_patter$folder_ud, "spatstat", "h"))

#### Define datasets
# NB: the detection data is used b/c assemble_acoustics() is called under the hood
datasets <- list(detections_by_unit = detections_by_unit, 
                 moorings = moorings,
                 archival_by_unit = archival_by_unit, 
                 behaviour_by_unit = behaviour_by_unit)

#### Estimate coordinates
# Implementation
iteration <- copy(iteration_patter)
lapply_estimate_coord_patter(iteration = iteration, datasets = datasets)
# Examine selected coords 
lapply_qplot_coord(iteration, 
                   "coord-fwd.qs",
                   extract_coord = function(s) s$states[sample.int(1000, size = .N, replace = TRUE), ])

#### Record
# Up to 55
# * BoundsError: attempt to access 250000-element Vector{Float64} at index [250001]
# - 50, 
# * [crop] cannot create dataset -> out of storage!

#### Estimate UDs
# We estimate UDs for iterations 1:3 with max number of particles
# * TO DO, filter iteration! 
# Implementation 
iteration[, file_coord := file.path(folder_coord, "coord-smo.qs")]
lapply_estimate_ud_spatstat(iteration = iteration, 
                            extract_coord = function(s) s$states,
                            cl = NULL, 
                            plot = FALSE)
# (optional) Examine selected UDs
lapply_qplot_ud(iteration, "spatstat", "h", "ud.tif")


###########################
###########################
#### Synthesis



#### End of code. 
###########################
###########################
