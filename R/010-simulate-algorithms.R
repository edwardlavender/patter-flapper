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
if (patter:::os_linux()) {
  stopifnot(!any(c("terra", "sf") %in% sort(loadedNamespaces())))
}

#### Load data 
if (!patter:::os_linux()) {
  map      <- terra::rast(here_data("spatial", "bathy.tif"))
  ud_grid  <- terra::rast(here_data("spatial", "ud-grid.tif"))
  ud_null  <- terra::rast(here_data("spatial", "ud-null.tif"))
  win      <- qs::qread(dv::here_data("spatial", "win.qs"))
} 
skateids <- qs::qread(here_data("input", "mefs", "skateids.qs"))
moorings <- qs::qread(here_data("input", "mefs", "moorings.qs")) 
pars     <- qs::qread(here_data("input", "pars.qs"))


###########################
###########################
#### Set up algorithms 

#### Julia set up
julia_connect()
set_seed()
set_model_move_components()

#### Define simulation timeline
n_path   <- 100L
timeline <- seq(as.POSIXct("2024-04-01 00:00:00", tz = "UTC"), 
                as.POSIXct("2024-04-30 23:58:00", tz = "UTC"), 
                by = "2 mins")

#### Define global datasets
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

# (optional) Cleanup
if (FALSE) {
  unlink(here_data("input", "simulation"), recursive = TRUE)
  dir.create(here_data("input", "simulation"))
  dirs.create(here_data("input", "simulation", seq_len(n_path)))
}
  
# Run simulation (~50 mins)
if (FALSE) {
  
  #### Simulate data on high-resolution grid
  map_5m <- terra::rast(here_data("spatial", "bathy-5m.tif"))
  set_map(map_5m)
  
  #### Define ModelMove structure to simulate paths
  # We simulate the path using the best-guess parameters
  model_move <- patter_ModelMove(pars$pmovement[1, ])
  
  #### Define ModelObs structures to simulate observations 
  # We simulate parameters using best-guess parameters
  # * These are hard-coded as 20, 20 for the depth observation model
  ModelObsAcousticLogisTruncPars <- 
    moorings |> 
    select(sensor_id = "receiver_id", "receiver_x", "receiver_y") |> 
    mutate(receiver_alpha = pars$pdetection$receiver_alpha[1], 
           receiver_beta = pars$pdetection$receiver_beta[1], 
           receiver_gamma = pars$pdetection$receiver_gamma[1]) |> 
    as.data.table()
  ModelObsDepthNormalTruncPars <- data.table(sensor_id = 1L, 
                                             depth_sigma = 20, 
                                             depth_deep_eps = 20)
  ModelObsPars <- list(ModelObsAcousticLogisTrunc = ModelObsAcousticLogisTruncPars, 
                       ModelObsDepthNormalTrunc = ModelObsDepthNormalTruncPars)
  
  # Iteratively simulate N movement paths 
  tic()
  path_id <- 1L
  while (path_id <= n_path) {
    
    print(path_id)

    #### Simulate initial 'tagging' location
    # We sample an initial location from the receiver array 
    # This helps to ensure we generate some detections for the COA/RSP/AC*PF algorithms
    xinit_bb  <- terra::ext(min(moorings$receiver_x), max(moorings$receiver_x), 
                            min(moorings$receiver_y), max(moorings$receiver_y))
    xinit_map <- terra::crop(map, xinit_bb)
    xinit_fwd <- terra::spatSample(xinit_map, size = 1L, xy = TRUE, na.rm = TRUE)
    xinit_fwd <- data.table(map_value = xinit_fwd$map_value,
                            x = xinit_fwd$x, 
                            y = xinit_fwd$y)
    if (model_move_is_crw()) {
      xinit_fwd[, angle := runif(.N) * 2 * pi]
    }
    
    #### Simulate behavioural states 
    # behaviour <- sample(c(1L, 2L), length(timeline), replace = TRUE)
    behaviour   <- simulate_behaviour(timeline)
    julia_assign("behaviour", behaviour)
    
    #### Define movement model
    state <- state_flapper
    update_model_move_components()
    
    #### Simulate movement path 
    coord_path <- sim_path_walk(.map = map, 
                                .timeline = timeline, 
                                .state = state, 
                                .xinit = xinit_fwd, 
                                .model_move = model_move, 
                                .n_path = 1L, 
                                .plot = FALSE)
    # points(moorings$receiver_x, moorings$receiver_y)
    # Record recapture location
    if (model_move_is_crw()) {
      xinit_bwd <- coord_path[.N, .(map_value, x, y, angle)]
    }
    
    #### Simulate observations
    yobs <- sim_observations(.timeline = timeline, 
                             .model_obs = ModelObsPars)
    # Check we have simulated detections
    ndet <- length(which(yobs$ModelObsAcousticLogisTrunc[[1]]$obs == 1L))
    if (ndet < 100L) {
      warn(paste(ndet, "detection(s) simulated."))
      next 
    }
    # Check depth time stamps are not duplicated
    if (any(duplicated(yobs$ModelObsDepthNormalTrunc[[1]]$timestamp))) {
      stop("Simulated depth time series contains duplicated time stamps.")
    }
    
    #### Map path UD 
    # * 13 s 500 pixels 
    stopifnot(spatstat.geom::spatstat.options("npixel") == 500)
    ud_path <- map_dens(.map        = ud_grid, 
                        .owin       = win,
                        .coord      = coord_path[, .(x, y)],
                        .discretise = TRUE,
                        sigma       = bw.h, 
                        .fterra     = TRUE, 
                        .plot       = FALSE)
    
    #### Save datasets
    here_out <- function(...) {
      here_data("input", "simulation", path_id, ...)
    }
    qs::qsave(xinit_fwd, here_out("xinit-fwd.qs"))
    qs::qsave(xinit_bwd, here_out("xinit-bwd.qs"))
    qs::qsave(behaviour, here_out("behaviour.qs"))
    qs::qsave(coord_path, here_out("coord.qs"))
    qs::qsave(yobs, here_out("yobs.qs"))
    terra::writeRaster(ud_path$ud, here_out("ud.tif"), overwrite = TRUE)
  
    #### Continue
    path_id <- path_id + 1L
    
  }
  toc()
  
  #### Visualise simulated paths (~33)
  mapfiles <- here_data("input", "simulation", seq_len(n_path), "ud.tif")
  mapfiles <- data.table(row = rep(1:10, each = 10),
                         column = rep(1:10, 10), 
                         mapfile = mapfiles)
  ggplot_maps(mapdt = mapfiles, 
              png_args = list(filename = here_fig("simulation", "map-paths.png"), 
                              height = 12, width = 15, units = "in", res = 600))
  
}

#### Load simulated datasets
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
  iteration <- copy(iteration_coa)
  if (FALSE) {
    
    #### Estimate coordinates (~2 s for three paths, ~83 s for 100 paths)
    iteration[, file_coord := file.path(folder_coord, "coord.qs")]
    lapply_estimate_coord_coa(iteration = iteration, datasets = datasets)
    # (optional) Examine selected coords
    lapply_qplot_coord(iteration, "coord.qs")
    
    #### Estimate UDs
    # Implementation (22 mins, 500 pixels, sigma = bw.h, cl = 10L)
    nrow(iteration)
    lapply_estimate_ud_spatstat(iteration = iteration, 
                                extract_coord = NULL,
                                cl = 1L, 
                                plot = FALSE)
    # (optional) Examine selected UDs
    lapply_qplot_ud(iteration, "spatstat", "h", "ud.tif")
    
  }
  
  #### Visualise ME
  # Compute ME 
  me <- pbapply::pbsapply(split(iteration, seq_row(iteration)), function(sim) {
    tryCatch(
      skill_me(.obs = terra::rast(here_data("input", "simulation", sim$unit_id, "ud.tif")), 
               .mod = terra::rast(file.path(sim$folder_ud, "spatstat", "h", "ud.tif"))), 
               error = function(e) NA)
  })
  # Compute ME for null model
  me_null <- pbapply::pbsapply(split(iteration, seq_row(iteration)), function(sim) {
    skill_me(.obs = terra::rast(here_data("input", "simulation", sim$unit_id, "ud.tif")), 
             .mod = ud_null)
  })
  # Visualise ME ~ delta_t
  iteration |>
    mutate(unit_id     = factor(unit_id), 
           delta_t     = factor(delta_t, levels = pcoa$delta_t), 
           me          = me,
    ) |> 
    filter(!is.na(me)) |>
    group_by(unit_id) |> 
    # Centre ME to facilitate plotting 
    mutate(me = me - mean(me)) |>
    ungroup() |> 
    as.data.table() |>
    ggplot() + 
    geom_point(aes(delta_t, me, colour = factor(unit_id), group = factor(unit_id))) + 
    geom_line(aes(delta_t, me, colour = factor(unit_id), group = factor(unit_id))) +
    geom_smooth(aes(as.integer(delta_t_int), me), method = "gam")  + 
    theme(legend.position = "none")
  
  # > Best guess:       2 days
  # > Restricted value: 1 day
  # > Flexible value:   3 days 
  
  ### Visualise maps
  # Define panel row (path) and column (parameter) labels
  cols <- c("NA", "2 days", "1 days", "3 days")
  cols <- factor(cols, levels = cols)
  rows <- c("Path", "COA[1]", "COA[2]", "COA[3]")
  rows <- factor(rows, levels = rows)
  panel <- data.table(row = rows, column = cols)
  # Define mapfiles for selected algorithm runs 
  selected_paths <- 1:3L
  mapfiles_alg <- 
    iteration |> 
    filter(unit_id %in% selected_paths) |>
    mutate(row = unit_id) |>
    filter(delta_t %in% panel$column) |> 
    mutate(column = panel$column[match(delta_t, panel$column)]) |>
    mutate(mapfile = file.path(folder_ud, "spatstat", "h", "ud.tif")) |>
    select(row, column, mapfile) |> 
    as.data.table()
  # Define mapfiles for simulated paths
  mapfiles_path <- 
    data.table(row = selected_paths, # n_path
               column = cols[1],
               mapfile = here_data("input", "simulation", selected_paths, "ud.tif")
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
  prsp <- data.table(parameter_id = 1:25L, 
                     er.ad = c(2.5,  # 250 * 0.01, 
                               12.5, # 250 * 0.05, # default
                               25.0,
                               seq(50, 1000, by = 50),
                               1250, 
                               1500
                     ))
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
  # unlink(iteration_rsp$folder_coord, recursive = TRUE)
  # unlink(iteration_rsp$folder_ud, recursive = TRUE)
  # list.files(iteration_rsp$folder_coord, recursive = TRUE)
  # list.files(iteration_rsp$folder_ud, recursive = TRUE)
  dirs.create(iteration_rsp$folder_coord)
  dirs.create(iteration_rsp$folder_ud)
  dirs.create(file.path(iteration_rsp$folder_ud, "spatstat", "h"))
  dirs.create(file.path(iteration_rsp$folder_ud, "dbbmm"))
  
  #### Define datasets
  datasets <- list(detections_by_unit = detections_by_unit, 
                   moorings = copy(moorings))
  
  #### Run algorithm
  iteration <- copy(iteration_rsp)
  if (FALSE) {
    
    #### Estimate coordinates (~89 mins)
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
    # Implementation (5 hr, 10 cl)
    lapply_estimate_ud_dbbmm(iteration = iteration, 
                             cl = 10L, 
                             plot = FALSE)
    # (optional) Examine selected UDs
    lapply_qplot_ud(iteration, "dbbmm", "ud.tif")
    # Load example maps with different er.ad values
    ud_500 <- terra::rast("data/output/simulation/1/rsp/ac/13/1/ud/dbbmm/ud.tif")
    ud_250 <- terra::rast("data/output/simulation/1/rsp/ac/8/1/ud/dbbmm/ud.tif")
    ud_750 <- terra::rast("data/output/simulation/1/rsp/ac/18/1/ud/dbbmm/ud.tif")
    # Visualise maps
    # * Visualise maps in sequence to visualise differences
    # * Maps are different as expected
    # * But er.ad parameter has very limited influence on results
    # * er.ad = 250 produces slightly more restricted map
    # * er.ad = 750 produces slightly more diffuse map
    # * ... in which low probability areas are spread out a bit more
    # * ... but differences are small. 
    # * There are also tiny changes in the size of the core region
    terra::plot(ud_250) 
    terra::plot(ud_500)
    terra::plot(ud_750) 
  }
  
  #### Visualise ME 
  me <- pbapply::pbsapply(split(iteration, seq_row(iteration)), function(sim) {
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
  iteration |>
    mutate(unit_id   = factor(unit_id), 
           er.ad     = factor(er.ad, levels = prsp$er.ad), 
           me        = me) |> 
    filter(!is.na(me)) |>
    group_by(unit_id) |> 
    # Centre ME to facilitate plotting 
    mutate(me = me - mean(me)) |>
    ungroup() |> 
    as.data.table() |>
    ggplot() + 
    geom_point(aes(er.ad, me, colour = factor(unit_id), group = factor(unit_id))) + 
    geom_line(aes(er.ad, me, colour = factor(unit_id), group = factor(unit_id))) +
    geom_smooth(aes(as.integer(er.ad), me), method = "gam") + 
    theme(legend.position = "none")
  
  # > Best guess:       ~500 
  # > Restricted value: ~250
  # > Flexible value:   ~750 
  
  ### Visualise maps
  # Define panel row (path) and column (parameter) labels
  cols <- c("NA", "500", "250", "750")
  cols <- factor(cols, levels = cols)
  rows <- c("Path", "RSP[1]", "RSP[2]", "RSP[3]")
  rows <- factor(rows, levels = rows)
  panel <- data.table(row = rows, column = cols)
  # Define mapfiles for selected algorithm runs 
  selected_paths <- 1:3L
  mapfiles_alg <- 
    iteration |>
    filter(unit_id %in% selected_paths) |>
    mutate(row = unit_id) |>
    filter(er.ad %in% panel$column) |> 
    mutate(column = panel$column[match(er.ad, panel$column)]) |>
    mutate(mapfile = file.path(folder_ud, "dbbmm", "ud.tif")) |>
    select(row, column, mapfile) |> 
    as.data.table()
  # Define mapfiles for simulated paths
  mapfiles_path <- 
    data.table(row = selected_paths, 
               column = cols[1],
               mapfile = here_data("input", "simulation", selected_paths, "ud.tif")
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

#### Define selected paths
selected_paths <- 1:3L

#### Define unitsets (unit_ids & algorithms)
unitsets <- 
  CJ(unit_id = selected_paths, 
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
# NB: only one value used for np now, this is set below
np <- NA_integer_
ni <- 1:3L
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
# Update np
nplookup <- data.table(dataset = c("ac", "dc", "acdc"), np = c(5000L, 5000L, 100000L))
iteration_patter[, np := nplookup$np[match(dataset, nplookup$dataset)]]
# sim$month_id is required by estimate_coord_patter() to define the timeline
iteration_patter[, month_id := "04-2024"]
# xinit folders
iteration_patter[, folder_xinit := file.path("data", "input", "simulation", unit_id)]
all(file.exists(file.path(iteration_patter$folder_xinit, "xinit-fwd.qs"))) |> stopifnot()
all(file.exists(file.path(iteration_patter$folder_xinit, "xinit-bwd.qs"))) |> stopifnot()
# Implement smoothing for all simulations
iteration_patter[, smooth := TRUE]

#### Build patter folders
dirs.create(iteration_patter$folder_coord)
dirs.create(iteration_patter$folder_ud)
dirs.create(file.path(iteration_patter$folder_ud, "spatstat", "h"))

#### Define datasets
# NB: the detection data is used b/c assemble_acoustics() is called under the hood
datasets <- list(detections_by_unit = copy(detections_by_unit[selected_paths]), 
                 moorings = moorings,
                 archival_by_unit = copy(archival_by_unit[selected_paths]), 
                 behaviour_by_unit = behaviour_by_unit[selected_paths])
# (essential) Process archival units for patter_ModelObs()
for (i in selected_paths) {
  datasets$archival_by_unit[[i]] <- 
    datasets$archival_by_unit[[i]] |> 
    select(timestamp, depth = obs) |> 
    # Delete parameter columns so they can be added by patter_ModelObs()
    mutate(depth_sigma = NULL, 
           depth_deep_eps = NULL) |> 
    as.data.table()
  
}
head(datasets$archival_by_unit[[1]])
head(datasets$archival_by_unit[[2]])

#### Initialise coordinate estimation 
# Define iterations 
iteration <- copy(iteration_patter)
# Set map
set_map(here_data("spatial", "bathy.tif"))
# Set batch & export vmap
# * We implement the algorithms in batches so that we export vmap once
batch     <- pars$pmovement$mobility[1]
iteration <- iteration[mobility == batch, ]
vmap      <- here_data("spatial", glue("vmap-{batch}.tif"))
set_vmap(.vmap = vmap)
rm(vmap)
# Check progress between loop restarts
if (FALSE) {
  it <- iteration[1:200, ]
  progress <- 
    lapply(split(it, seq_row(it)), function(d) {
      # d <- it[1, ]
      qs::qread(file.path(d$folder_coord, "data-fwd.qs"))
    }) |> 
    rbindlist(fill = TRUE)
}
# (optional) Test convergence
if (FALSE) {
  # Test convergence for selected algorithm
  # * AC: 5000 particles: success for 1:3
  # * DC: 5000 particles: success for 1:3
  # * ACDC: 
  # - 100,000: 1 (fail); 2 (success); 3 (TO DO), ... 
  # - 200,000: 1 (fail)
  # - 250,000: 1 (success: 17 mins); 
  
  # iteration_trial <- iteration[dataset == "acdc" & sensitivity == "best" & iter == 1L, ]
  iteration_trial <- iteration_patter[iter == 1L, ]
  iteration_trial <- arrange(iteration_trial, dataset, sensitivity)
  iteration_trial[dataset == "acdc", np := 250000] 
  iteration_trial[, smooth := FALSE]
  iteration_trial <- iteration_trial[dataset == "ac", ]
  # iteration_trial <- iteration_trial[1, ]
  # debug(estimate_coord_patter)
  nrow(iteration_trial)
  lapply_estimate_coord_patter(iteration = iteration_trial, 
                               datasets = datasets,
                               trial = TRUE)
  qs::qread(file.path(iteration$folder_coord[1], "data-fwd.qs"))
  # Compare output
  if (!patter:::os_linux()) {
    here_data("input", "simulation", "1", "ud.tif") |> terra::rast() |> terra::plot()
    map_pou(.map = terra::rast(here_data("input", "spatial", "map.tif")),
            .coord = file.path(iteration$folder_coord[1], "coord-smo.qs"))
  }

}
# Implementation
gc()
nrow(iteration)
lapply_estimate_coord_patter(iteration = iteration, datasets = datasets)
# Examine selected coords 
if (patter:::os_linux()) {
  stop("Exit server at this point (for convenience).")
}
lapply_qplot_coord(iteration, 
                   "coord-fwd.qs",
                   extract_coord = function(s) s$states[sample.int(1000, size = .N, replace = TRUE), ])

#### Estimate UDs
# We estimate UDs for iterations 1:3 with max number of particles
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