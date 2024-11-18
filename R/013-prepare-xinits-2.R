###########################
###########################
#### prepare-xinits-2.R

#### Aims
# (1) Prepare initial locations for the filter (part 2)
# * This code samples initial locations for the filter:
# * It is important that initial locations are as precise as possible
# * Simulations show that this facilitates convergence
# * We identified the nearest known location (recapture event, detection)
# * We simulated 30 realisations of the prior from that location in the intervening time
# * We compute the maximum radius from the last known position (min: 2000 m, plus buffer if long gap)
# * We sample initial locations within this region 
# * We assume at the starting point that the individual is ±20 m of high-resolution seabed data
# * These restrictions form an 'improper prior'
# * They are more precise than strictly implied by the prior but really help the filter. 
# * We use the same initial locations for all analyses (ACDC, DC, AC) to facilitate convergence. 

#### Prerequisites
# (1) Prepare datasets (recaps, acoustics, archival)


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
# try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
dv::src()
library(Tools4ETS)

#### Load data
# Load map (high-resolution)
map        <- terra::rast(here_data("spatial", "bathy-5m.tif"))
unitsets   <- qs::qread(here_data("input", "unitsets.qs"))
moorings   <- qs::qread(here_data("input", "mefs", "moorings.qs"))
detections <- qs::qread(here_data("input", "mefs", "acoustics.qs"))
archival   <- qs::qread(here_data("input", "mefs", "archival.qs"))
recaps     <- qs::qread(here_data("input", "mefs", "recaps.qs"))
pars       <- qs::qread(here_data("input", "pars.qs"))


###########################
###########################
#### Data preparation 

#### Julia connect
julia_connect()
set_seed()
julia_command('Rasters.checkmem!(false);')
set_map(map)

#### Aggregate map
# The aggregated map is used to improve plotting speed only 
background <- terra::aggregate(map, fact = 20, fun = "mean", na.rm = TRUE)

#### Prepare unitsets
# Define unique unit_ids and corresponding start dates
# For each unit_id (individual, month) we will use the best possible starting values
unitsets <- 
  unitsets |>
  group_by(unit_id) |> 
  slice(1L) |> 
  mutate(timestamp = paste0("01-", month_id), 
         timestamp = as.Date(timestamp, "%d-%m-%Y"), 
         timestamp = as.POSIXct(paste(timestamp, "00:00:00"), tz = "UTC")) |> 
  select(unit_id, individual_id, month_id, timestamp) |> 
  as.data.frame()

#### Prepare recaps
recaps    <- as.data.frame(recaps)

#### Prepare detections
detections <- 
  detections |> 
  left_join(moorings, by = "receiver_id") |> 
  select(individual_id, timestamp, x = receiver_x, y = receiver_y) |> 
  # Add buffer column
  # * (Detection time stamps are known accurately, unlike some recapture time stamps)
  mutate(buffer = 0) |>
  as.data.frame()

#### Prepare archival 
# Round time steps 
# This is necessary for one of the checks below
archival[, timestamp := lubridate::round_date(timestamp, "2 mins")]
archival  <- as.data.frame(archival)
# stopifnot(archival$timestamp[2376076] == as.POSIXct("2016-04-01 00:00:00", tz = "UTC"))


###########################
###########################
#### Simulate starting locations

#### Clean up
directions <- c("forward", "backward")
if (FALSE) {
  unlink(here_data("input", "xinit"), recursive = TRUE)
}
sapply(c("forward", "backward"), function(direction) {
  folders <- here_data("input", "xinit", rho, unitsets$individual_id, unitsets$month_id)
  dirs.create(folders)
}) |> invisible()

#### Define movement model (RW or CRW)
set_model_move_components()
model_move <- patter_ModelMove(pars$pmovement[1, ])

#### Simulate starting locations (~7.5 mins per direction)
lapply(directions, function(direction) {

  tic()
  for (i in 1:nrow(unitsets)) {
    
    # direction = "forward"
    # i = 1
    cat("\n")
    print(paste0(rep("-", 50), collapse = ""))
    print(paste(i, "/", nrow(unitsets)))
    
    #### Define individual data
    # For direction = "forward", sim$timestamp contains the start time
    # For direction = "backward", we update sim$timestamp to the end of the month
    sim <- unitsets[i, ]
    if (direction == "backward") {
      sim$timestamp <- sim$timestamp + months(1) - 120
      print(sim$timestamp)
    }
    
    #### Identify nearest observations 
    # Positions 
    precap <- match_ts_nearest_by_key(sim, recaps, "individual_id", "timestamp")
    pdet   <- match_ts_nearest_by_key(sim, detections, "individual_id", "timestamp")
    parc   <- match_ts_nearest_by_key(sim, archival, "individual_id", "timestamp")
    # Corresponding time stamps 
    trecap <- recaps$timestamp[precap]
    tdet   <- detections$timestamp[pdet]
    tarc   <- archival$timestamp[parc]
    
    #### Identify nearest known coordinate 
    # Compute time difference before fixed-location observation 
    gap <- difftime(sim$timestamp, c(trecap, tdet), units = "mins")  |> abs() |> as.numeric()
    # Identify corresponding coordinate from recapture event or acoustic observation
    if (which.min(gap) == 1) {
      gap   <- gap[1]
      coord <- recaps[precap, ]
    } else {
      gap   <- gap[2]
      coord <- detections[pdet, ]
    }
    if (coord$buffer > 0) {
      warning("Using coordinates with uncertain timing (buffer > 0). Check 'nearest observation' assignment!", immediate. = TRUE)
      print(trecap); print(tdet); print(tarc)
    }
    # Define xinit for simulation 
    xinit <- 
      coord |> 
      mutate(map_value = terra::extract(map, cbind(x, y))) |> 
      select(map_value, x, y) |> 
      as.data.table()
    # Include angle for CRW
    if (model_move_is_crw()) {
      xinit[, angle := runif(.N) * 2 * pi]
    }
    
    #### Compute reasonable maximum radius of individual from point according to prior 
    radius <- NULL
    if (gap <= 2) {
      radius <- 2000 
    }
    if (is.null(radius)) {
      # Define timeline 
      timeline    <- seq(as.POSIXct("2016-01-01 00:00:00", tz = "UTC"), by = "2 mins", length.out = ceiling(gap / 2))
      # Iteratively simulate movements from point & compute distance travelled from point 
      radius <- 
        pbapply::pbsapply(1:30, function(j) {
          
          # Simulate behaviour & set movement model 
          state       <- state_flapper
          behaviour   <- simulate_behaviour(timeline)
          behaviour   <- simulate_behaviour(timeline)
          julia_assign("behaviour", behaviour)
          update_model_move_components()
          
          # Simulate path 
          path <- 
            sim_path_walk(.map = map, 
                          .timeline = timeline, 
                          .state = state, 
                          .xinit = xinit, 
                          .model_move = model_move, 
                          .n_path = 1L, 
                          .plot = FALSE) |> 
            suppressWarnings()
          
          # Compute travelled radius from point  
          radius <- terra::distance(cbind(xinit$x, xinit$y), cbind(path$x, path$y), lonlat = FALSE)[1, ]
          max(radius)
          
        }) |> max() 
      
      # (optional) Add buffer to radius for more prolonged gaps
      if (gap > 60 * 12) {
        radius <- radius + 5000
      }
    }
    radius <- max(c(2000, radius))
    radius <- radius + coord$buffer
    stopifnot(length(radius) == 1L)
    
    #### Define region within which to sample initial locations 
    # Define possible locations of individual 
    container <- 
      cbind(xinit$x, xinit$y) |>
      terra::vect(crs = terra::crs(map)) |>
      terra::buffer(width = radius)
    # terra::plot(map)
    # terra::lines(container)
    # Define origin SpatRaster
    # * NB: mask = TRUE has a bug: we implement masking in a separate step
    origin <- terra::crop(map, container, mask = FALSE)
    origin <- terra::mask(origin, container)
    # terra::plot(origin)
    # Restrict origin by depth 
    if (sim$timestamp == tarc) {
      depth   <- archival$depth[parc] 
      deep    <- origin + 20
      shallow <- origin - 20
      msk <- deep >= depth & shallow <= depth
      # terra::plot(msk)
      origin <- terra::mask(origin, msk, maskvalues = FALSE)
      # terra::plot(origin)
    } else {
      warn("No depth observation for row!")
    }
    
    #### Sample initial locations 
    # Sample locations
    # NB: as.data.frame(origin, xy = TRUE, na.rm = TRUE may exhaust vector memory)
    xinit <- terra::spatSample(origin, size = 1e6L, replace = FALSE, 
                               na.rm = TRUE, values = TRUE, xy = TRUE)
    if (nrow(xinit) == 0L) {
      terra::plot(map)
      terra::lines(container)
      abort("No initial locations for row!")
    }
    # Process data.table
    xinit <-
      xinit |> 
      select(map_value, x, y) |> 
      as.data.table()
    # Add turning angle, if applicable
    if (model_move_is_crw()) {
      xinit[, angle := runif(.N) * 2 * pi]
    }
    
    #### Write xinit to file
    if (direction == "forward") {
      xinit.qs <- "xinit-fwd.qs"
    } else if (direction == "backward") {
      xinit.qs <- "xinit-bwd.qs"
    }
    qs::qsave(xinit, here_data("input", "xinit", rho, sim$individual_id, sim$month_id, xinit.qs))
    
    invisible(NULL)
    
  }
  toc()
  NULL
  
})

# NB: recapture events with uncertain timing are never used.

###########################
###########################
#### Visualise xinits

#### Iterate over directions & plot xinits (~1.6 mins)
tic()
lapply(c("fwd", "bwd"), function(direction) {
  
  # List files 
  xinits <- list.files(here_data("input", "xinit", rho), recursive = TRUE, full.names = TRUE)
  xinits <- gtools::mixedsort(xinits)
  xinits <- xinits[stringr::str_detect(xinits, direction)]
  
  # Make figure
  png(here_fig("xinit", glue("xinits-{rho}-{direction}.png")), 
      height = 12, width = 12, units = "in", res = 600)
  pp <- par(mfrow = par_mf(length(xinits)))
  cl_lapply(xinits, function(xinit) {
    # Read initial locations 
    unit_id <- xinit |> tools::file_path_sans_ext() |> basename()
    xinit   <- qs::qread(xinit)
    # Plot low-resolution background (for speed) plus xinits
    terra::plot(background, 
                main = unit_id, font = 2,
                legend = FALSE)
    points(xinit$x, xinit$y, pch = ".")
    NULL
  })
  dev.off()
  nothing()
}) |> invisible()
toc()


#### End of code. 
###########################
###########################