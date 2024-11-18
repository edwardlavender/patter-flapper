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
  win      <- qs::qread(here_data("spatial", "win.qs"))
  mpa      <- qreadvect(here_data("spatial", "mpa.qs"))
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
# set_ModelObsCaptureContainer()

#### Define simulation timeline
n_path         <- 100L
selected_paths <- 1:3L
timeline       <- seq(as.POSIXct("2024-04-01 00:00:00", tz = "UTC"), 
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

#### Flag moorings in MPA
# The moorings_in_mpa dataset is used to compute detection days
if (!patter:::os_linux()) {
  
  # Visualise receivers in MPA
  terra::plot(map)
  terra::plot(mpa, "id")
  points(moorings$receiver_x, moorings$receiver_y, pch = 4)
  
  # Identify moorings in MPA 
  moorings_in_mpa <- 
    moorings |>
    mutate(open = terra::extract(mpa, cbind(moorings$receiver_x, moorings$receiver_y))$open, 
           open = as.character(open)) |>
    filter(!is.na(open)) |>
    as.data.table()
  
  # All moorings within MPA are in closed zones
  table(moorings_in_mpa$open)
  
}


###########################
###########################
#### Simulate paths and observations

# (optional) Cleanup
if (FALSE) {
  unlink(here_data("input", "simulation"), recursive = TRUE)
  list.files(here_data("input", "simulation"), recursive = TRUE)
  dir.create(here_data("input", "simulation"))
  dirs.create(here_data("input", "simulation", seq_len(n_path)))
}
  
# Run simulation (~45 mins)
if (FALSE) {
  
  #### Simulate data on selected grid
  # 1) Use aggregated grid and correct observation model
  # - This tells us how well algorithms perform if we get parameters right
  # 2) Use high-resolution grid for simulation & aggregated grid for modelling
  # - This indicates whether our understanding 
  # - ... of the underlying process + aggregation effects is correct
  # - But this requires many more particles to achieve convergence.
  # map_5m <- terra::rast(here_data("spatial", "bathy-5m.tif"))
  # set_map(map_5m)
  set_map(map)
  
  #### Define ModelMove structure to simulate paths
  # We simulate the path using the best-guess parameters
  model_move <- patter_ModelMove(pars$pmovement[1, ])
  
  #### Define ModelObs structures to simulate observations 
  # We simulate parameters using best-guess parameters
  # * These may be hard-coded for the high-resolution grid (20, 20)
  ModelObsAcousticLogisTruncPars <- 
    moorings |> 
    select(sensor_id = "receiver_id", "receiver_x", "receiver_y") |> 
    mutate(receiver_alpha = pars$pdetection$receiver_alpha[1], 
           receiver_beta = pars$pdetection$receiver_beta[1], 
           receiver_gamma = pars$pdetection$receiver_gamma[1]) |> 
    as.data.table()
  # ModelObsDepthNormalTruncPars <- data.table(sensor_id = 1L, 
  #                                            depth_sigma = 20, 
  #                                            depth_deep_eps = 20)
  ModelObsDepthNormalTruncPars <- data.table(sensor_id = 1L, 
                                             depth_sigma = pars$pdepth$depth_sigma[1], 
                                             depth_deep_eps = pars$pdepth$depth_deep_eps[1])
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
    coord_path <- sim_path_walk(.map        = map, 
                                .timeline   = timeline, 
                                .state      = state, 
                                .xinit      = xinit_fwd, 
                                .model_move = model_move, 
                                .n_path     = 1L, 
                                .plot       = FALSE)
    xinit_bwd <- coord_path[.N, .(map_value, x, y)]
    # points(moorings$receiver_x, moorings$receiver_y)
    
    #### Record capture/recapture locations for filter
    # We assume the capture/recapture locations are known
    # Angles are unknown
    xinit_fwd <- lapply(1:1e5L, function(i) {
      xinit_fwd
    }) |> rbindlist()
    if (model_move_is_crw()) {
      xinit_fwd[, angle := runif(.N) * 2 * pi]
    }
    xinit_bwd <- lapply(1:1e5L, function(i) {
      xinit_bwd
    }) |> rbindlist()
    if (model_move_is_crw()) {
      xinit_bwd[, angle := runif(.N) * 2 * pi]
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
              png_args = list(filename = here_fig("simulation", "map-paths-full.png"), 
                              height = 12, width = 15, units = "in", res = 600))
  
  #### Visualise selected paths
  mapfiles <- here_data("input", "simulation", selected_paths, "ud.tif")
  mapfiles <- data.table(row = c(1, 2, 3),
                         column = c(1, 1, 1), 
                         mapfile = mapfiles)
  ggplot_maps(mapdt = mapfiles, 
              png_args = list(filename = here_fig("simulation", "map-paths.png"), 
                              height = 12, width = 15, units = "in", res = 600))
  
  #### Estimate 'true' residency (~10 s)
  iteration_res <- data.table(unit_id = seq_len(n_path), 
                              algorithm = "Path", 
                              sensitivity = "Best", 
                              iter = 1L,
                              file = here_data("input", "simulation", seq_len(n_path), "coord.qs"))
  residency <- 
    lapply_estimate_residency_coord(files = iteration_res$file,
                                    cl = 10L)
  residency <- 
    left_join(residency, iteration_res, by = "file") |>
    mutate(file = NULL) |> 
    select(unit_id, algorithm, sensitivity, iter, zone, time) |> 
    arrange(unit_id, algorithm, sensitivity, iter, zone) |>
    as.data.table()
  qs::qsave(residency, here_data("output", "simulation-summary", "residency-path.qs"))
  
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
# Recapture containers
# xinit_by_unit <- lapply(seq_len(n_path), function(path) {
#   list(
#     forward  = qs::qread(here_data("input", "simulation", path, "xinit-fwd.qs")),
#     backward = qs::qread(here_data("input", "simulation", path, "xinit-bwd.qs"))
#   )
# })


###########################
###########################
#### Null model

if (FALSE) {
  
  #### Estimate map skill (ME) 
  me <- cl_lapply(seq_len(n_path), .cl = 10L, .fun = function(unit_id) {
    skill_me(.obs = terra::rast(here_data("input", "simulation", unit_id, "ud.tif")), 
             .mod = ud_null)
  }) |> unlist() |> as.numeric()
  me <- data.table(unit_id = seq_len(n_path), 
                   algorithm = "Null", 
                   sensitivity = "Best", 
                   iter = 1L, 
                   me = me)
  qs::qsave(me, here_data("output", "simulation-summary", "me-null.qs"))
  
  #### Estimate residency under null model
  study_area <- 
    terra::expanse(ud_grid)[, 2]
  residency_null_closed <- 
    sum(terra::expanse(mpa)[which(mpa$open == "closed")]) / study_area
  residency_null_open <- 
    sum(terra::expanse(mpa)[which(mpa$open == "open")]) / study_area
  residency_null_total <- 
    sum(terra::expanse(mpa)) / study_area
  residency_null <- 
    lapply(seq_len(n_path), function(unit_id) {
      data.table(unit_id = unit_id, 
                 algorithm = "Null", 
                 sensitivity = "Best", 
                 iter = 1L, 
                 zone = c("closed", "open", "total"), 
                 time = c(residency_null_closed, residency_null_open, residency_null_total))
    }) |> rbindlist()
  qs::qsave(residency_null, here_data("output", "simulation-summary", "residency-null.qs"))
  
}


###########################
###########################
#### Detection days

#### Estimate residency from detection days
if (FALSE & !patter:::os_linux()) {
  
  residency <- lapply(seq_len(n_path), function(path) {
    
    # Total number of days in month 
    ndays <- as.integer(lubridate::days_in_month(timeline[1]))
    
    # Compute detection days for receivers in MPA
    dds_total <- 
      detections_by_unit[[path]] |> 
      filter(receiver_id %in% moorings_in_mpa$receiver_id) |> 
      mutate(day = lubridate::day(timestamp)) |> 
      summarise(time = length(unique(day)) / ndays) |>
      mutate(unit_id = path, 
             algorithm = "DD", 
             sensitivity = "Best",
             iter = 1L,
             zone    = "total") |> 
      select(unit_id, algorithm, sensitivity, iter, zone, time) |> 
      arrange(unit_id, algorithm, sensitivity, iter, zone) |>
      as.data.table()
    
    # Detection days in closed areas are identical 
    # * (All receivers were in closed areas)
    # * For the figures, we also record DDs in closed areas
    dds_closed <- copy(dds_total)
    dds_closed[, zone := "closed"]
    
    rbind(dds_total, dds_closed)
    
  }) |> rbindlist()
  
  qs::qsave(residency, 
            here_data("output", "simulation-summary", "residency-detection-days.qs"))
  
}

###########################
###########################
#### COA algorithm

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
  if (FALSE) {
    unlink(iteration_coa$folder_coord, recursive = TRUE)
    unlink(iteration_coa$folder_ud, recursive = TRUE)
    list.files(iteration_coa$folder_coord, recursive = TRUE)
    list.files(iteration_coa$folder_ud, recursive = TRUE)
  }
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
                                cl = 10L, 
                                plot = FALSE)
    # (optional) Examine selected UDs
    lapply_qplot_ud(iteration, "spatstat", "h", "ud.tif")
    
  }
  
  #### Select delta_t
  # Compute ME
  me <- cl_lapply(split(iteration, seq_row(iteration)), .cl = 10L, .fun = function(sim) {
    tryCatch(
      skill_me(.obs = terra::rast(here_data("input", "simulation", sim$unit_id, "ud.tif")), 
               .mod = terra::rast(file.path(sim$folder_ud, "spatstat", "h", "ud.tif"))), 
               error = function(e) NA)
  }) |> unlist() |> as.numeric()
  # Compute ME for null model
  me_null <- cl_lapply(split(iteration, seq_row(iteration)), .cl = 10L, .fun = function(sim) {
    skill_me(.obs = terra::rast(here_data("input", "simulation", sim$unit_id, "ud.tif")), 
             .mod = ud_null)
  }) |> unlist() |> as.numeric()
  # Visualise ME ~ delta_t
  png(here_fig("simulation", "parameterisation-coa.png"), 
      height = 4, width = 5, units = "in", res = 600)
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
    geom_point(aes(delta_t, me, colour = factor(unit_id), 
                   group = factor(unit_id)), alpha = 0.5, size = 0.8) + 
    geom_line(aes(delta_t, me, colour = factor(unit_id), 
                  group = factor(unit_id)), alpha = 0.4, linewidth = 0.5) +
    geom_smooth(aes(as.integer(delta_t), me), method = "gam", 
                colour = "black", linewidth = 1.2, fill = "dimgrey")  + 
    theme(legend.position = "none") + 
    scale_y_continuous(limits = c(-4e-6, 6e-6), expand = c(0, 0), labels = prettyGraphics::sci_notation) + 
    xlab(expression(Delta * T)) + 
    ylab(expression("Relative mean error (" * italic(RME) * ")")) + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.position = "none")
  dev.off()
  # Select delta_t
  # > Best guess:       2 days
  # > Restricted value: 1 day
  # > Flexible value:   3 days 
  selected_delta_t <- c("2 days", "1 day", "3 days")
  stopifnot(all(selected_delta_t %in% iteration$delta_t))
  
  #### Visualise maps for selected_paths
  # Define panel row (path) and column (parameter) labels
  cols <- c("Path", selected_delta_t)
  cols <- factor(cols, levels = cols)
  rows <- c("Path", "COA[1]", "COA[2]", "COA[3]")
  rows <- factor(rows, levels = rows)
  panel <- data.table(row = rows, column = cols)
  # Define mapfiles for selected algorithm runs 
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
  
  #### Estimate ME (across all paths)
  me <- 
    iteration |> 
    mutate(algorithm = "COA", 
           sensitivity = factor(delta_t, levels = selected_delta_t, labels = c("Best", "AC(-)", "AC(+)")), 
           me = me) |>
    filter(delta_t %in% selected_delta_t) |> 
    select(unit_id, algorithm, sensitivity, iter, me) |> 
    arrange(unit_id, algorithm, sensitivity, iter) |>
    as.data.table()
  # Check convergence for 'best' simulations (100 %)
  table(is.na(me$me))
  me |>
    filter(sensitivity == "Best") |> 
    summarise(convergence_pr = length(which(!is.na(me))) / n())
  qs::qsave(me, here_data("output", "simulation-summary", "me-coa.qs"))
  
  #### Estimate residency (across all paths)
  # Define dataset
  iteration_res <- 
    iteration |>
    # filter(unit_id %in% selected_paths) |> 
    filter(delta_t %in% selected_delta_t) |> 
    mutate(file = file.path(folder_ud, "spatstat", "h", "ud.tif")) |> 
    as.data.table()
  # Compute residency 
  residency <- lapply_estimate_residency_ud(files = iteration_res$file)
  # Write output
  residency <- 
    left_join(iteration_res, residency, by = "file") |> 
    mutate(algorithm = "COA", 
           sensitivity = factor(delta_t, levels = selected_delta_t, labels = c("Best", "AC(-)", "AC(+)"))) |>
    select(unit_id, algorithm, sensitivity, iter, zone, time) |> 
    arrange(unit_id, algorithm, sensitivity, iter, zone) |>
    as.data.table()
  qs::qsave(residency, here_data("output", "simulation-summary", "residency-coa.qs"))

}


###########################
###########################
#### RSP algorithm

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
  if (FALSE) {
    unlink(iteration_rsp$folder_coord, recursive = TRUE)
    unlink(iteration_rsp$folder_ud, recursive = TRUE)
    list.files(iteration_rsp$folder_coord, recursive = TRUE)
    list.files(iteration_rsp$folder_ud, recursive = TRUE)
  }
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
    
    #### Estimate coordinates (~1.5 hours)
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
  
  #### Check errors & computation time
  # Larger er.ad reduces mean error (no clear optimum) (see below)
  # But there are costs:
  # 1) Errors become more frequent, even with an grid the size of Scotland
  # 2) Computation time increases
  # Therefore, in the selection of the 'best' RSP values, we have to be pragmatic:
  # * Select a sufficiently large value to reduce ME
  # * But not so large as to impose errors & increased computation time 
  # * ... with only minimal continued improvement in ME
  rspinfo <- cl_lapply(split(iteration, seq_row(iteration)), .fun = function(sim) {
    file <- file.path(sim$folder_ud, "dbbmm", "data.qs")
    if (file.exists(file)) {
      out <- qs::qread(file)
      out[success == FALSE, time := NA_real_]
      out[, er.ad := sim$er.ad]
    } else {
      out <- NULL
    }
    return(out)
  }) |> rbindlist()
  # We observe errors with er.ad >= 950 m: 
  unique(rspinfo$error)
  min(rspinfo$er.ad[!is.na(rspinfo$error)])
  # Errors are apparent at larger values:
  yjit <- runif(nrow(rspinfo), 0, 0.1)
  plot(rspinfo$er.ad, as.numeric(!is.na(rspinfo$error)) + yjit)
  # Computation time grows, on average, with er.ad:
  plot(rspinfo$er.ad, rspinfo$time)
  
  #### Select er.ad
  # Compute ME
  me <- cl_lapply(split(iteration, seq_row(iteration)), .cl = 10L, .fun = function(sim) {
    tryCatch(   
      skill_me(.obs = terra::rast(here_data("input", "simulation", sim$unit_id, "ud.tif")), 
               .mod = terra::rast(file.path(sim$folder_ud, "dbbmm", "ud.tif"))), 
      error = function(e) NA)
  }) |> unlist() |> as.numeric()
  me_null <- cl_lapply(split(iteration, seq_row(iteration)), .cl = 10L, .fun = function(sim) {
      skill_me(.obs = terra::rast(here_data("input", "simulation", sim$unit_id, "ud.tif")), 
               .mod = ud_null)
  }) |> unlist() |> as.numeric()
  # Visualise ME ~ er.ad
  png(here_fig("simulation", "parameterisation-rsp.png"), 
      height = 4, width = 5, units = "in", res = 600)
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
    geom_point(aes(er.ad, me, colour = factor(unit_id), 
                   group = factor(unit_id)), alpha = 0.5, size = 0.8) + 
    geom_line(aes(er.ad, me, colour = factor(unit_id), 
                  group = factor(unit_id)), alpha = 0.4, linewidth = 0.5) +
    geom_smooth(aes(as.integer(er.ad), me), method = "gam", 
                colour = "black", linewidth = 1.2, fill = "dimgrey")  + 
    theme(legend.position = "none") + 
    scale_y_continuous(limits = c(-4e-6, 6e-6), expand = c(0, 0), labels = prettyGraphics::sci_notation) + 
    xlab(expression(er.ad)) + 
    ylab(expression("Relative mean error (" * italic(RME) * ")")) + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title.x = element_text(family = "mono"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.position = "none")
  dev.off()
  # Select er.ad
  # > Best guess:       ~500 
  # > Restricted value: ~250
  # > Flexible value:   ~750 
  selected_er.ad <- c(500, 250, 750)
  stopifnot(all(selected_er.ad %in% iteration$er.ad))
  
  #### Visualise maps for selected_paths
  # Define panel row (path) and column (parameter) labels
  cols <- c("Path", selected_er.ad)
  cols <- factor(cols, levels = cols)
  rows <- c("Path", "RSP[1]", "RSP[2]", "RSP[3]")
  rows <- factor(rows, levels = rows)
  panel <- data.table(row = rows, column = cols)
  # Define mapfiles for selected algorithm runs 
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
  
  #### Estimate ME (across all paths)
  me <- 
    iteration |> 
    mutate(algorithm = "RSP", 
           sensitivity = factor(er.ad, levels = selected_er.ad, labels = c("Best", "AC(-)", "AC(+)")), 
           me = me) |>
    filter(er.ad %in% selected_er.ad) |>
    select(unit_id, algorithm, sensitivity, iter, me) |> 
    arrange(unit_id, algorithm, sensitivity, iter) |>
    as.data.table()
  # Check convergence for 'best' simulations (100 %)
  table(is.na(me$me))
  me |>
    filter(sensitivity == "Best") |> 
    summarise(convergence_pr = length(which(!is.na(me))) / n())
  qs::qsave(me, here_data("output", "simulation-summary", "me-rsp.qs"))
  
  #### Estimate residency (across all paths)
  # Define dataset
  iteration_res <- 
    iteration |>
    # filter(unit_id %in% selected_paths) |> 
    filter(er.ad %in% selected_er.ad) |> 
    mutate(file = file.path(folder_ud, "dbbmm", "ud.tif")) |> 
    as.data.table()
  # Compute residency 
  residency <- lapply_estimate_residency_ud(files = iteration_res$file)
  # Write output
  residency <- 
    left_join(iteration_res, residency, by = "file") |> 
    mutate(algorithm = "RSP", 
           sensitivity = factor(er.ad, levels = selected_er.ad, labels = c("Best", "AC(-)", "AC(+)"))) |>
    select(unit_id, algorithm, sensitivity, iter, zone, time) |> 
    arrange(unit_id, algorithm, sensitivity, iter, zone) |>
    as.data.table()
  qs::qsave(residency, here_data("output", "simulation-summary", "residency-rsp.qs"))
  
}


###########################
###########################
#### Patter algorithms

#### Define unitsets (unit_ids & algorithms)
unitsets <- 
  CJ(unit_id = seq_len(n_path), # selected_paths 
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
  
  # d <- unitsets[1, ]
  
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
  
  # Combine data with all parameters
  d <- cbind(d, p)
  
  # Select relevant rows
  # * For selected_paths, we use all parameters
  # * Otherwise, we just use 1 iteration & the best parameters
  if (!(d$unit_id[1] %in% selected_paths)) {
     d <- d[iter == 1L & sensitivity == "best", ]
  }
  
  d
  
}) |> 
  rbindlist() |> 
  mutate(folder = file.path("data", "output", "simulation", unit_id, "patter"),
         folder_coord = file.path(folder, dataset, parameter_id, iter, "coord"), 
         folder_ud = file.path(folder, dataset, parameter_id, iter, "ud")) |>
  distinct() |> 
  arrange(unit_id, dataset, parameter_id) |>
  mutate(index = row_number()) |>
  dplyr::select("index", 
                "unit_id",
                "dataset", 
                "parameter_id", "iter", "sensitivity", 
                "k1", "k2", "theta1", "theta2", "mobility", 
                "receiver_alpha", "receiver_beta",
                "receiver_gamma", "depth_sigma", "depth_deep_eps",
                "np",
                "folder_coord", "folder_ud") |>
  distinct() |>
  as.data.table()
# Check rows (153 -> 444)
nrow(iteration_patter)
# View(iteration_patter)

#### Add supporting columns
# Update np
nplookup <- data.table(dataset = c("ac", "dc", "acdc"), np = c(5000L, 5000L, 10000L))
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
if (FALSE) {
  unlink(iteration_patter$folder_coord, recursive = TRUE)
  unlink(iteration_patter$folder_ud, recursive = TRUE)
  list.files(iteration_patter$folder_coord, recursive = TRUE)
  list.files(iteration_patter$folder_ud, recursive = TRUE)
}
dirs.create(iteration_patter$folder_coord)
dirs.create(iteration_patter$folder_ud)
dirs.create(file.path(iteration_patter$folder_ud, "spatstat", "h"))

#### Define datasets
# NB: the detection data is used b/c assemble_acoustics() is called under the hood
datasets <- list(detections_by_unit = copy(detections_by_unit), 
                 moorings = moorings,
                 archival_by_unit = copy(archival_by_unit), 
                 behaviour_by_unit = behaviour_by_unit)
# (Essential) Process archival units for patter_ModelObs()
for (i in seq_len(length(datasets$archival_by_unit))) {
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

#### Run algorithm
iteration <- copy(iteration_patter)
if (FALSE) {
  
  #### Initialise coordinate estimation 
  # Set map
  set_map(here_data("spatial", "bathy.tif"))
  # Set batch & export vmap
  # * We implement the algorithms in batches so that we export vmap once
  batch     <- pars$pmovement$mobility[1]
  iteration <- iteration[mobility == batch, ]
  vmap      <- here_data("spatial", glue("vmap-{batch}.tif"))
  set_vmap(.vmap = vmap)
  rm(vmap)
  # (optional) Select rows/dataset
  # * batch 1: 390 rows (tmux sessions 1:8)
  # * batch 2: 27 rows  (tmux session 9)
  # * batch 3: 27 rows  (tmux session 10)
  rows <- seq_row(iteration)
  if (batch == pars$pmovement$mobility[1]) {
    rows <- 1:50      # 1
    # rows <- 51:100  # 2 
    # rows <- 101:150 # 3
    # rows <- 151:200 # 4
    # rows <- 201:250 # 5
    # rows <- 251:300 # 6
    # rows <- 301:350 # 7
    # rows <- 351:390 # 8 
  }
  iteration <- iteration[rows, ]
  
  #### (optional) Check progress between loop restarts
  if (FALSE) {
    it <- iteration[1:200, ]
    progress <- 
      lapply(split(it, seq_row(it)), function(d) {
        # d <- it[1, ]
        qs::qread(file.path(d$folder_coord, "data-fwd.qs"))
      }) |> 
      rbindlist(fill = TRUE)
  }
  
  #### (optional) Test convergence
  if (FALSE) {
    
    # Test convergence for selected algorithm (best)
    # * AC    : 5000 particles    : 100 % success 
    # * DC    : 5000 particles    : 92 % success
    # * ACDC  : 10000 particles   : 98 % success
    
    # Approximate wall time per run
    # * MacBook: 
    # - 2 mins x 2 (5,000 particles) for filter
    # - 15-20 mins (500 particles) for smoother 
    
    # Select dataset of interest
    iteration_trial <- copy(iteration_patter)
    iteration_trial <- iteration_trial[iter == 1L, ]
    iteration_trial <- iteration_trial[sensitivity == "best", ]
    # iteration_trial <- iteration_trial[unit_id == 1, ]
    iteration_trial <- arrange(iteration_trial, dataset, sensitivity)
    iteration_trial <- iteration_trial[dataset == "ac", ]
    # iteration_trial <- iteration_trial[1, ]
    iteration_trial[, np := 5000L] 
    iteration_trial[, smooth := FALSE]
    # debug(estimate_coord_patter)
    nrow(iteration_trial)
    lapply_estimate_coord_patter(iteration = iteration_trial, 
                                 datasets = datasets,
                                 trial = TRUE, 
                                 log.folder = here_data("output", "log", "simulation", "trials"),
                                 log.txt = "ac.txt")
    qs::qread(file.path(iteration_trial$folder_coord[1], "data-fwd.qs"))
    
    # Compare output
    if (!patter:::os_linux()) {
      here_data("input", "simulation", "1", "ud.tif") |> terra::rast() |> terra::plot()
      map_pou(.map = terra::rast(here_data("input", "spatial", "map.tif")),
              .coord = file.path(iteration$folder_coord[1], "coord-smo.qs"))
    }
    
  }
  
  #### Estimate coords
  gc()
  nrow(iteration)
  head(iteration)
  lapply_estimate_coord_patter(iteration = iteration, 
                               datasets = datasets, 
                               trial = FALSE, 
                               log.folder = here_data("output", "log", "simulation"), 
                               log.txt = glue("{batch}-{min(rows)}.txt"))
  
  # Examine selected coords (slow)
  # if (patter:::os_linux()) {
  #   stop("Exit server at this point (for convenience).")
  # }
  # lapply_qplot_coord(iteration,
  #                    "coord-smo.qs",
  #                    extract_coord = function(s) s$states[sample.int(1000, size = .N, replace = TRUE), ])
  
  #### Estimate UDs
  if (!patter:::os_linux()) {
    
    #### Clean up
    # Clean up folders
    if (FALSE) {
      unlink(file.path(iteration$folder_ud, "pou"), recursive = TRUE)
      unlink(file.path(iteration$folder_ud, "spatstat", "h"), recursive = TRUE)
    }
    # Rebuild
    dirs.create(file.path(iteration$folder_ud, "pou"))
    dirs.create(file.path(iteration$folder_ud, "spatstat", "h"))
    # Define input coordinates
    iteration[, file_coord := file.path(folder_coord, "coord-smo.qs")]
    
    #### Estimate UDs via POU (~7 mins)
    lapply_estimate_ud_pou(iteration = iteration,
                           extract_coord = function(s) s$states,
                           cl = 12L,
                           plot = FALSE)
    # (optional) Examine selected UDs
    lapply_qplot_ud(iteration, "pou", "ud.tif")
    
    #### Estimate UDs via spatstat (~23 mins)
    lapply_estimate_ud_spatstat(iteration = iteration,
                                extract_coord = function(s) s$states,
                                cl = 12L,
                                plot = FALSE)
    # (optional) Examine selected UDs
    lapply_qplot_ud(iteration, "spatstat", "h", "ud.tif")
  }

}

# Checks
iteration |> 
  filter(unit_id == 1L & dataset == "dc" & iter == 1L) |>
  select(parameter_id, sensitivity, folder_coord, folder_ud) |> 
  as.data.table()

list.files("data/output/simulation/1/patter/dc/2/1/coord")
list.files("data/output/simulation/1/patter/dc/2/1/ud", recursive = TRUE)


###########################
###########################
#### Patter analysis 

###########################
#### Tidy iteration 

if (FALSE) {
  
  #### Tidy iteration
  # Copy
  iteration <- copy(iteration_patter)
  # Define grouping variables with tidy labels
  iteration[, unit_id := factor(unit_id)]
  iteration[, dataset := factor(dataset, levels = c("ac", "dc", "acdc"), labels = c("AC", "DC", "ACDC"))]
  # Revise 'sensitivity' coding
  iteration[sensitivity == "movement" & mobility == pars$pmovement$mobility[2], sensitivity := "movement-under"]
  iteration[sensitivity == "movement" & mobility == pars$pmovement$mobility[3], sensitivity := "movement-over"]
  iteration[sensitivity == "ac" & receiver_alpha == pars$pdetection$receiver_alpha[2], sensitivity := "ac-under"]
  iteration[sensitivity == "ac" & receiver_alpha == pars$pdetection$receiver_alpha[3], sensitivity := "ac-over"]
  iteration[sensitivity == "dc" & depth_sigma == pars$pdepth$depth_sigma[2], sensitivity := "dc-under"]
  iteration[sensitivity == "dc" & depth_sigma == pars$pdepth$depth_sigma[3], sensitivity := "dc-over"]
  iteration[, sensitivity := factor(sensitivity, 
                                    levels = c("best", 
                                               "movement-under", "movement-over", 
                                               "ac-under", "ac-over",
                                               "dc-under", "dc-over"), 
                                    labels = c("Best", 
                                               "Move (-)", "Move (+)", 
                                               "AC(-)", "AC(+)", 
                                               "DC(-)", "DC(+)"))] 
  
  algorithms     <- c("Null", "DD", "COA", "RSP", "AC", "DC", "ACDC")
  col_all        <- rainbow(length(algorithms), alpha = 0.3)
  names(col_all) <- algorithms
  col_heuristic  <- col_all[1:4]
  col_patter     <- col_all[5:7]
  
}


###########################
#### Patter error check 

if (FALSE) {
  
  # Check for errors on forward filter 
  sapply(split(iteration, seq_row(iteration)), function(d) {
    qs::qread(file.path(d$folder_coord, "data-fwd.qs"))$error
  }) |> unlist() |> unique()
  
  # Check for errors on backward filter 
  sapply(split(iteration, seq_row(iteration)), function(d) {
    file <- file.path(d$folder_coord, "data-bwd.qs")
    if (file.exists(file)) {
      qs::qread(file)$error
    }
  }) |> unlist() |> unique()
  
  # Check for errors on smoother
  sapply(split(iteration, seq_row(iteration)), function(d) {
    file <- file.path(d$folder_coord, "data-smo.qs")
    if (file.exists(file)) {
      qs::qread(file)$error
    }
  }) |> unlist() |> unique()
  
}


###########################
#### Patter convergence  

if (FALSE) {
  
  #### Define convergence
  iteration[, convergence := sapply(split(iteration, seq_row(iteration)), function(d) {
    file.exists(file.path(d$folder_coord, "coord-smo.qs"))
  })]
  # Review overall convergence
  table(iteration$convergence)
  # Review convergence for ACDC
  iteration[dataset == "ACDC" & sensitivity == "Best", .(unit_id, iter, convergence)]
  
  #### Visualise convergence (best)
  # This is a barplot of the proportion of simulations that converged
  # (If we had multiple iterations for each path, it could be a boxplot)
  png(here_fig("simulation", "convergence-best.png"), 
      height = 4, width = 3, units = "in", res = 600)
  iteration |> 
    filter(sensitivity == "Best" & iter == 1L) |> 
    group_by(dataset) |> 
    summarise(convergence = length(which(convergence)) / n()) |>
    ungroup() |>
    as.data.table() |> 
    ggplot(aes(dataset, convergence, fill = dataset)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = col_patter) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
    labs(
      x = "Algorithm",
      y = "Pr(convergence)",
      fill = "Algorithm"
    ) + 
    theme_bw() + 
    theme(panel.spacing.y = unit(1.5, "lines"), 
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 10)),
          panel.grid = element_blank())
  dev.off()
  
  #### Visualise convergence (sensitivity)
  # Plot Pr(convergence) ~ algorithm with panels for best/under/over estimation of parameters
  # This plot is only useful if convergence is not 100 %
  png(here_fig("simulation", "convergence-sensitivity.png"), 
      height = 4.5, width = 9, units = "in", res = 600)
  iteration |>
    filter(unit_id %in% selected_paths) |>
    group_by(unit_id, dataset, sensitivity) |> 
    summarise(convergence = length(which(convergence)) / n()) |>
    as.data.table() |> 
    ggplot(aes(x = dataset, ymin = 0, ymax = convergence, colour = dataset)) + 
    geom_linerange(linewidth = 1.2) +
    scale_colour_manual(values = col_patter) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
    labs(
      x = "Algorithm",
      y = "Pr(convergence)",
      colour = "Algorithm"
    ) + 
    facet_grid(unit_id ~ sensitivity, 
               scales = "free_x") +
    theme_bw() + 
    theme(panel.spacing.y = unit(1.5, "lines"), 
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 10)),
          panel.grid.minor.y = element_blank(), 
          panel.grid.major.y = element_blank())
  dev.off()
  
}


###########################
#### Patter maps

if (FALSE) {
  
  #### For each selected path, visualise UDs for performance
  # Define mapfiles for simulated path
  mapfiles_path <- 
    data.table(row = selected_paths, 
               column = "Path",
               mapfile = here_data("input", "simulation", selected_paths, "ud.tif")
    )
  # Select iterations  
  mapfiles_alg <- 
    iteration |> 
    filter(unit_id %in% selected_paths) |>
    filter(sensitivity == "Best") |> 
    filter(iter == 1L) |> 
    mutate(row = unit_id, 
           column = dataset, 
           mapfile = file.path(folder_ud, "spatstat", "h", "ud.tif")) |> 
    select(row, column, mapfile) |>
    as.data.table()
  # Collect mapfiles (row, column, mapfile)
  mapfiles <- 
    rbind(mapfiles_alg, mapfiles_path) |> 
    arrange(row, column) |> 
    mutate(row = factor(row, levels = selected_paths), 
           column = factor(column, levels = c("Path", "AC", "DC", "ACDC"))) |>
    as.data.table()
  # Make maps
  ggplot_maps(mapdt = mapfiles, 
              png_args = list(filename = here_fig("simulation", "map-patter-spatstat.png"), 
                              height = 5, width = 10, units = "in", res = 600))
  
  #### For each selected path, visualise UDs for performance & sensitivity (~15 s)
  # (Manually simulated path maps with patter UD maps outside of R)
  cl_lapply(selected_paths, function(path) {
    
    #### Define iteration dataset for mapping
    # path <- selected_paths[1]
    iteration_map <- iteration[iter == 1L, ]
    iteration_map <- iteration_map[unit_id == path, ]
    iteration_map[, key := paste(dataset, sensitivity)]
    if (FALSE) {
      # Visual check of simulated and reconstructed patter for path 1 
      pp <- par(mfrow = c(1, 2))
      terra::plot(terra::rast("data/input/simulation/1/ud.tif"))
      terra::plot(terra::rast("data/output/simulation/1/patter/ac/1/1/ud/spatstat/h/ud.tif"))
      par(pp)
    }

    #### Define file paths to ud.tif (POU or spatstat)
    # Use POU
    # udtype <- "pou"
    # iteration_map[, mapfile := file.path(iteration_map$folder_ud, "pou", "ud.tif")]
    # Use spatstat
    udtype <- "spatstat"
    iteration_map[, mapfile := file.path(iteration_map$folder_ud, "spatstat", "h", "ud.tif")]
    
    #### Define mapfiles data.table for mapping
    mapfiles <- CJ(row = levels(iteration_map$dataset), 
                   column = levels(iteration_map$sensitivity), 
                   mapfile = NA)
    mapfiles[, row := factor(row, levels = levels(iteration_map$dataset))]
    mapfiles[, column := factor(column, levels = levels(iteration_map$sensitivity))]
    mapfiles[, key := paste(row, column)]
    mapfiles[, mapfile := iteration_map$mapfile[match(key, iteration_map$key)]]
    # Drop NA panels
    mapfiles <- mapfiles[!(key %in% c("AC DC(-)", "AC DC(+)", "DC AC(-)", "DC AC(+)")), ]
    
    #### Make maps
    ggplot_maps(mapdt = mapfiles, 
                png_args = list(filename = here_fig("simulation", glue("map-patter-{udtype}-{path}.png")), 
                                height = 3, width = 7, units = "in", res = 600))
    
  })
 
  #### For one selected path, for each algorithm, visualise UDs for repeatability
  cl_lapply(c("AC", "DC", "ACDC"), function(alg) {
    
    #### Define iteration dataset for mapping
    # alg <- "AC"
    iteration_map <- iteration[unit_id == 1L, ]
    iteration_map <- iteration_map[dataset == alg, ]
    iteration_map[, iter := factor(iter, levels = sort(unique(iter)), labels = paste("Iter", sort(unique(iter))))]
    iteration_map[, key := paste(iter, sensitivity)]
    
    #### Define file paths to ud.tif (POU or spatstat)
    # Use POU
    # udtype <- "pou"
    # iteration_map[, mapfile := file.path(iteration_map$folder_ud, "pou", "ud.tif")]
    # Use spatstat
    udtype <- "spatstat"
    iteration_map[, mapfile := file.path(iteration_map$folder_ud, "spatstat", "h", "ud.tif")]
    
    #### Define mapfiles data.table for mapping
    mapfiles <- CJ(row = unique(iteration_map$iter), 
                   column = levels(iteration_map$sensitivity), 
                   mapfile = NA)
    mapfiles[, row := factor(row, levels = levels(iteration_map$iter))]
    mapfiles[, column := factor(column, levels = levels(iteration_map$sensitivity))]
    mapfiles[, key := paste(row, column)]
    mapfiles[, mapfile := iteration_map$mapfile[match(key, iteration_map$key)]]
    # Drop NA panels
    if (alg == "AC") {
      mapfiles <- mapfiles[!(column %in% c("DC(-)", "DC(+)")), ]
    }
    if (alg == "DC") {
      mapfiles <- mapfiles[!(column %in% c("AC(-)", "AC(+)")), ]
    }
    
    #### Make maps
    ggplot_maps(mapdt = mapfiles, 
                png_args = list(filename = here_fig("simulation", glue("map-patter-{udtype}-{alg}.png")), 
                                height = 3, width = 7,units = "in", res = 600))
    
  })
   
}


###########################
#### Patter ME

if (FALSE) {
  
  # Compute ME
  me <- cl_lapply(split(iteration, seq_row(iteration)), .cl = 10L, .fun = function(sim) {
    # ud.tif <- file.path(sim$folder_ud, "pou", "ud.tif")
    ud.tif <- file.path(sim$folder_ud, "spatstat", "h", "ud.tif")
    tryCatch(   
      skill_me(.obs = terra::rast(here_data("input", "simulation", sim$unit_id, "ud.tif")), 
               .mod = terra::rast(ud.tif)), 
      error = function(e) NA)
  }) |> unlist() |> as.numeric()
  
  # Collate ME
  me <-
    iteration |>
    mutate(algorithm = dataset,
           me = me) |>
    select(unit_id, algorithm, sensitivity, iter, me) |>
    arrange(unit_id, algorithm, sensitivity, iter) |>
    as.data.table()
  
  # qs::qsave(me, here_data("output", "simulation-summary", "me-patter-pou.qs"))
  qs::qsave(me, here_data("output", "simulation-summary", "me-patter-spatstat.qs"))

}


###########################
#### Patter residency 

if (FALSE) {
  
  # Compute residency (~2 min, 10 cl)
  iteration_res <- copy(iteration)
  iteration_res[, file := file.path(folder_coord, "coord-smo.qs")]
  residency <- lapply_estimate_residency_coord(files = iteration_res$file, 
                                               extract_coord = function(s) s$states, 
                                               cl = 10L)
  # Write output
  residency <- 
    left_join(iteration_res, residency, by = "file") |> 
    mutate(algorithm = dataset) |>
    select(unit_id, algorithm, sensitivity, iter, zone, time) |> 
    arrange(unit_id, algorithm, sensitivity, iter, zone) |>
    as.data.table()
  qs::qsave(residency, here_data("output", "simulation-summary", "residency-patter.qs"))
  
}


###########################
###########################
#### Synthesis

# TO DO
# Unify colour scales

###########################
#### Maps synthesis

#### Combine ME datasets
skill <- 
  rbindlist(
    list(
      qs::qread(here_data("output", "simulation-summary", "me-null.qs")),
      qs::qread(here_data("output", "simulation-summary", "me-coa.qs")),
      qs::qread(here_data("output", "simulation-summary", "me-rsp.qs")),
      # qs::qread(here_data("output", "simulation-summary", "me-patter-pou.qs"))
      # Use spatstat for improved representation of patter patterns
      qs::qread(here_data("output", "simulation-summary", "me-patter-spatstat.qs"))
    )
  )

#### Visualise overall map skill (ME ~ algorithm)
png(here_fig("simulation", "me.png"), 
    height = 3, width = 6, units = "in", res = 600)
skill |> 
  filter(sensitivity == "Best") |>
  ggplot() + 
  geom_violin(aes(algorithm, me, fill = algorithm), 
              linewidth = 0.25, 
              scale = "count") +
  scale_y_continuous(
    limits = c(0, 1.7e-5),
    expand = c(0, 0),
    labels = prettyGraphics::sci_notation) + 
  scale_fill_manual(values = col_all) + 
  xlab("Algorithm") + 
  ylab(expression("Mean absolute error (" * italic(ME) * ")")) + 
  labs(fill = "Algorithm") +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) 
dev.off()

#### Visualise map skill sensitivity for selected_paths
# For Null, COA, RSP  : 1 point per panel
# For AC, DC, ACDC    : 3 points per panel (three iterations)
png(here_fig("simulation", "me-sensitivity.png"), 
    height = 5, width = 12, units = "in", res = 600)
skill |> 
  filter(unit_id %in% selected_paths) |> 
  as.data.table() |>
  ggplot() + 
  geom_jitter(aes(algorithm, me, fill = algorithm), 
              width = 0.2,
             shape = 21) + 
  scale_fill_manual(values = col_all) + 
  scale_y_continuous(labels = prettyGraphics::sci_notation) + 
  xlab("Algorithm") + 
  ylab(expression("Mean absolute error (" * italic(ME) * ")")) + 
  labs(fill = "Algorithm") +
  facet_grid(unit_id ~ sensitivity, scales = "free_x") + 
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)), 
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


###########################
#### Residency synthesis 

if (FALSE) {
  
  #### Combine residency datasets
  residency <- 
    rbindlist(
      list(
        qs::qread(here_data("output", "simulation-summary", "residency-path.qs")),
        qs::qread(here_data("output", "simulation-summary", "residency-null.qs")),
        qs::qread(here_data("output", "simulation-summary", "residency-detection-days.qs")),
        qs::qread(here_data("output", "simulation-summary", "residency-coa.qs")),
        qs::qread(here_data("output", "simulation-summary", "residency-rsp.qs")),
        qs::qread(here_data("output", "simulation-summary", "residency-patter.qs"))
      )
    )
  
  #### Examine structure
  head(residency)
  unique(residency$algorithm)
  unique(residency$sensitivity)
  residency <- 
    residency |>
    mutate(zone = factor(zone, 
                         levels = c("open", "closed", "total"), 
                         labels = c("Open", "Closed", "Protected")), 
           sensitivity = factor(sensitivity, 
                                levels = c("Best", 
                                           "Move (-)", "Move (+)",
                                           "AC(-)", "AC(+)", 
                                           "DC(-)", "DC(+)"))) |>
    arrange(unit_id, algorithm, sensitivity, iter, zone) |> 
    as.data.table()

  #### Visualise overall residency skill (by zone)
  # Define residency estimates for best-case simulations 
  residency_skill <- 
    residency |> 
    filter(sensitivity == "Best") |> 
    # filter(zone == "total") |> 
    as.data.table()
  # Define true residency
  residency_skill_path <- 
    residency_skill |> 
    filter(algorithm == "Path") |>
    select(unit_id, zone, truth = time) |>
    as.data.table()
  # Define estimated residency 
  residency_skill_alg <- 
    residency_skill |>
    filter(algorithm != "Path") |> 
    select(unit_id, zone, algorithm, time) |>
    as.data.table()
  # Define error (estimated time - truth)
  # * <0: underestimation of residency
  # * >0: overestimation of residency  
  residency_skill <- 
    merge(residency_skill_alg, residency_skill_path, by = c("unit_id", "zone")) |> 
    mutate(error = time - truth) |> 
    as.data.table()
  # Visualise skill
  png(here_fig("simulation", "residency-mpa-skill.png"), 
      height = 3, width = 10, units = "in", res = 600)
  residency_skill |>
    ggplot() + 
    geom_violin(aes(algorithm, error, fill = algorithm), 
                linewidth = 0.25,
                scale = "count") + 
    geom_hline(yintercept = 0, linetype = 3) + 
    scale_fill_manual(values = col_all) + 
    scale_y_continuous(expand = c(0, 0), limits = c(-1, 1)) + 
    xlab("Algorithm") + 
    ylab(expression("Residency error (" * italic(RE) * ")")) + 
    labs(fill = "Algorithm") +
    facet_wrap(~zone, scales = "free_x") +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank(), 
          panel.grid.major.y = element_blank(), 
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 10))) 
  dev.off()
  
  # > Null model underestimates residency 
  # > DD underestimates residency b/c individuals linger around receivers (in our simulations) 
  # > COAs + RSPs overestimate residency b/c they don't permit movement away from receivers
  # > AC, DC, ACDC (TO DO)
  
  #### Visualise residency sensitivity for selected_paths
  png(here_fig("simulation", "residency-mpa-sensitivity.png"), 
      height = 6, width = 8, units = "in", res = 600)
  # Get true residency
  true_res <- 
    residency |>
    filter(unit_id %in% selected_paths) |>
    filter(algorithm == "Path", sensitivity == "Best") |>
    select(unit_id, zone, true_time = time) |> 
    as.data.table()
  # Pair true residency with 'best' estimate for arrows 
  arrow_res <- 
    true_res |> 
    merge(  
      residency |>
        filter(unit_id %in% selected_paths, sensitivity == "Best") |>
        filter(algorithm != "Path") |> 
        mutate(algorithm = factor(algorithm, levels = algorithms)) |> 
        select(unit_id, zone, algorithm, best_time = time) |> 
        as.data.table(), 
      by = c("unit_id", "zone"))
  # Make ggplot 
  residency |>
    filter(unit_id %in% selected_paths) |> 
    filter(algorithm != "Path") |>
    mutate(algorithm = factor(algorithm, levels = algorithms)) |>
    as.data.table() |>
    ggplot(aes(algorithm, time)) + 
    # Draw lines from the 'best' estimates to the truth (for visual clarity)
    geom_segment(
      data = arrow_res,
      aes(x = algorithm, xend = algorithm, y = true_time, yend = best_time),
      arrow = arrow(length = unit(0.15, "cm")),
      colour = "dimgrey",
      size = 0.5) + 
    # Add estimates with jitter
    geom_jitter(
      width = 0.25,
      shape = 21, size = 1.5, stroke = 0.5, 
      colour = "black", alpha = 0.6,
      aes(algorithm, time, fill = sensitivity)
    ) + 
    # Add best estimates as a horizontal line for visual clarity
    geom_hline(data = true_res,
               aes(yintercept = true_time),
               colour = "black") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
    xlab("Algorithm") + 
    ylab("Residency") + 
    labs(fill = "Parameterisation") + 
    facet_grid(unit_id ~ zone, scales = "free_x") + 
    # facet_grid(unit_id ~ sensitivity, scales = "free_x") + 
    theme_bw() + 
    theme(panel.spacing.y = unit(1.5, "lines"), 
          panel.grid.minor.y = element_blank(), 
          panel.grid.major.y = element_blank(), 
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 10)), 
          axis.text.x = element_text(angle = 45, hjust = 1))
  dev.off()
  
  # > Results
  # > TO DO
  
  # Summary statistics for performance simulations 
  residency_skill |> 
    group_by(zone, algorithm) |> 
    summarise(median = median(error, na.rm = TRUE), 
              se     = sd(error, na.rm = TRUE),
              lwr    = quantile(error, 0.025, na.rm = TRUE),
              upr    = quantile(error, 0.975, na.rm = TRUE)) |> 
    ungroup() |>
    prettyGraphics::tidy_numbers(digits = 2) |>
    as.data.table()
  # > COA median error: 10 %, 18 % uncertainty
  # > RSP median error: 9 %, 17 % uncertainty
  # > AC, DC and ACDC median error: 
  # > 1-2 % (5-10 fold improvement in median error) 
  # > reduced variation (1.5-3.6 fold reduction in uncertainty)
  
  # Summary statistics for sensitivity simulations
  # > Do the above results still hold if we consider sensitivity? 
  # > I.e., does the variation among paths envelope the variation among sensitivity trials? 
  
  
}


#### End of code. 
###########################
###########################