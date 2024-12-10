###########################
###########################
#### prepare-xinits-1.R

#### Aims
# (1) Prepare initial locations for the filter (part 1)
# * We use recapture events & acoustic detections to define xinit
# * This code prepares the recaptures dataset

#### Prerequisites
# (1) Download of skatespotter dataset (Lacie_Share)
# (2) Copy archival_raw.rds dataset    (Lacie_Share)


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
# try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
dv::src()
library(plotly)

#### Load data
map       <- terra::rast(here_data("spatial", "bathy.tif"))
skateids  <- qs::qread(here_data("input", "mefs", "skateids.qs"))
recaps    <- readRDS(here_data_raw("movement", "recaptures_processed.rds"))
archival  <- readRDS(here_data_raw("movement", "archival_raw.rds"))
moorings  <- qs::qread(here_data("input", "mefs", "moorings.qs")) 
iteration <- qs::qread(here_data("input", "iteration", "patter.qs"))
pars      <- qs::qread(here_data("input", "pars.qs"))


###########################
###########################
#### Assemble recaps

# Identify relevant skateids 
# (We focused on individuals with acoustic & archival data)
skateids <- skateids[skateids$individual_id %in% iteration$individual_id, ]

# Identify relevant recaptures 
recaps <- as.data.table(recaps)
recaps <- recaps[source %in% c("skateids", "rcs_database"), ]
recaps <- recaps[dst_id %in% skateids$dst_id, ]
recaps[, individual_id := skateids$individual_id[match(recaps$dst_id, skateids$dst_id)]]
recaps <- recaps[individual_id %in% skateids$individual_id, ]

# Isolate datasets by source
recaps_skateids <- recaps[source == "skateids", ]
recaps_skatespotter <- recaps[source == "rcs_database", ]

# For recaps_skateids, add deployment locations 
recaps_skateids[, long := skateids$long_tag_capture[match(individual_id, skateids$individual_id)]]
recaps_skateids[, lat := skateids$lat_tag_capture[match(individual_id, skateids$individual_id)]]

# Convert to UTM
xy <- cbind(recaps_skateids$long, recaps_skateids$lat) |> 
  terra::vect(crs = "WGS84") |> 
  terra::project(terra::crs(map)) |> 
  terra::crds()

# Tidy 
recaps_skateids <-
  recaps_skateids |>
  mutate(timestamp = as.POSIXct(date, tz = "UTC"),
         location = NA_character_,
         x = xy[, 1], 
         y = xy[, 2]) |> 
  select("individual_id", "timestamp", "x", "y") |>
  as.data.table()

# (optional) Check & assign locations
if (FALSE) {
  for (i in 1:nrow(recaps_skateids)) {
    terra::plot(map, main = paste(i))
    points(xy[i, 1], xy[i, 2], col = "red", lwd = 10)
    readline("Press [enter] to continue...")
  }
}
recaps_skateids[, location := "Kerrera"]
recaps_skateids[7, location := "Insh"]
if (FALSE) {
  for (i in 1:nrow(recaps_skateids)) {
    terra::plot(map, main = paste(i, recaps_skateids$location[i]))
    points(xy[i, 1], xy[i, 2], col = "red", lwd = 10)
    readline("Press [enter] to continue...")
  }
}

# For recaps_skatespotter, write xlsx & manually look up locations 
recaps_skatespotter |> 
  mutate(pit = skateids$pit_tag_id[match(individual_id, skateids$individual_id)],
         sex = skateids$sex[match(individual_id, skateids$individual_id)],
         timestamp = as.POSIXct(date, tz = "UTC")) |> 
  select(individual_id, dst_id, acoustic_id, pit, timestamp, long, lat) |>
  as.data.table() |> 
  xlsx::write.xlsx(here_data_raw("movement", "recaps-skatespotter.xlsx"))

# Read updated dataset
recaps_skatespotter <- 
  here_data_raw("movement", "recaps-skatespotter-updated.xlsx") |>
  xlsx::read.xlsx(1) |> 
  filter(!is.na(location)) |> 
  mutate(timestamp = as.POSIXct(timestamp, tz = "UTC")) |> 
  as.data.table()

# Convert coordinates to UTM 
xy <- cbind(recaps_skatespotter$long, recaps_skatespotter$lat) |> 
  terra::vect(crs = "WGS84") |> 
  terra::project(terra::crs(map)) |> 
  terra::geom(df = TRUE) |> 
  select("x", "y") 
recaps_skatespotter[, x := xy[, 1]]
recaps_skatespotter[, y := xy[, 2]]

# Validate locations: ok.
if (FALSE) {
  for (i in 1:nrow(recaps_skatespotter)) {
    terra::plot(map, main = paste(i, recaps_skatespotter$location[i]))
    points(xy[i, 1], xy[i, 2], col = "red", lwd = 10)
    readline("Press [enter] to continue...")
  }
}

# Set missing locations 
recaps_skateids
recaps_skatespotter[is.na(x) & location == "Kerrera", x := 708886.3]
recaps_skatespotter[is.na(y) & location == "Kerrera", y := 6254404]

# Collate recapture information
recaps <- 
  rbind(
    recaps_skateids |> select("individual_id", "timestamp", "x", "y"),
    recaps_skatespotter |> select("individual_id", "timestamp", "x", "y")
  ) |>
  arrange(individual_id, timestamp) |>
  as.data.table()


###########################
###########################
#### Refine time stamps

#### Method
# The above time stamps are Dates
# We record the dataset in Excel 
# We visualise the raw archival data around each time stamp
# We identify the exact time of time stamps

#### Write recaps spreadsheet to temporary file
recaps |> 
  xlsx::write.xlsx(here_data_raw("movement", "recaps.xlsx"))

#### Process raw archival time series for examination 
arc <- 
  archival |>
  mutate(individual_id = skateids$individual_id[match(dst_id, skateids$dst_id)], 
         timestamp = date_time,
         depth = abs(depth) * -1,
         temp = NULL) |> 
  filter(individual_id %in% recaps$individual_id) |>
  select(individual_id, timestamp, depth) |>
  arrange(individual_id, timestamp) |>
  as.data.frame()
stopifnot(nrow(arc) > 0L)

#### Manually examine depth time series
# And update the spreadsheet by examining depth time series around recapture events
if (FALSE) {
  
  recaps[, index := 1:.N]
  lapply(split(recaps, seq_row(recaps)), function(d) {
    
    # Define archival data for the recapture event
    # * Select individual
    # * Zoom in to window around event
    # d <- recaps[8, ]
    message(d$index)
    arc_for_recap <- 
      arc |> 
      filter(individual_id == d$individual_id) |> 
      filter(timestamp >= d$timestamp - 24 * 60 * 60 & timestamp <= d$timestamp + 24 * 60 * 60) |> 
      as.data.table()
    
    # Visualise time series, if available
    if (nrow(arc_for_recap) > 0L) {
      print(arc_for_recap[1, ])
      p <- 
        arc_for_recap |> 
        ggplot(aes(timestamp, depth)) + 
        geom_line() + 
        geom_point(size = 0.5, colour = "royalblue") + 
        ylim(c(-250, 0)) +
        labs(title = d$index)
      theme_bw()
      p <- ggplotly(p)  
      print(p)
    } else {
      # If the depth time series ended before the recapture event, we'll assume 14:00:00
      # * This is half-way between 08:00 and 20:00:00 (realistic bounds on recapture times)
      # * This is relevant for individual's 13 and 32
      cat("> Use $timestamp 12:00:00...")
    }
    
    readline(prompt = "Press [Enter] to continue...")
    NULL
    
  }) |> invisible()
  # prettyGraphics::vis_ts(arc)
  
}


###########################
###########################
#### Set buffer

#### Method
# We do not know the exact time of some recapture events
# We assume 14:00:00 (see above)
# This is ± 6 hours from the true recapture time 
# We simulate how far an individual can move in 8 hours
# And add this as an additional buffer column
# (However, recapture events with buffers are not required in prepare-xinits.R)

if (FALSE) {
  
  #### Initialise simulation 
  
  # Connect to Julia
  julia_connect()
  set_seed()
  set_model_move_components()
  set_map(map)
  
  # Define timeline
  # * We assume the set recapture time is within 8 hours of the true time 
  timeline <- seq(as.POSIXct("2024-04-01 08:00:00", tz = "UTC"), 
                  as.POSIXct("2024-04-01 14:00:00", tz = "UTC"), 
                  by = "2 mins")
  
  # Set movement model
  model_move <- patter_ModelMove(pars$pmovement[1, ])
  
  #### Simulate movement paths
  paths <- cl_lapply(1:30L, function(i) {
    
    
    #### Define movement model
    state       <- state_flapper
    behaviour   <- simulate_behaviour(timeline)
    julia_assign("behaviour", behaviour)
    update_model_move_components()
    
    #### Simulate movement path 
    # Simulate initial location in array 
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
    # Simulate coordinates 
    coord <- sim_path_walk(.map        = map, 
                           .timeline   = timeline, 
                           .state      = state, 
                           .xinit      = xinit_fwd, 
                           .model_move = model_move, 
                           .n_path     = 1L, 
                           .plot       = FALSE)
    # Update path_id 
    coord[, path_id := i]
    coord
  })
  
  #### Compute distance travelled in 8 hours
  # This value is set as the buffer in recaps-updated.xlsx
  paths |> 
    rbindlist() |>
    group_by(path_id) |>
    mutate(dist = terra::distance(matrix(c(x[1], y[1]), ncol = 2), cbind(x, y), lonlat = FALSE)[1, ]) |> 
    ungroup() |> 
    summarise(max = max(dist)) |>
    as.data.table()
  # max
  # <num>
  #   1: 5605.219
  
}

###########################
###########################
#### Collate recaptures


#### Collate recaps-updated.xlsx
recaps <- 
  here_data_raw("movement", "recaps-updated.xlsx") |>
  xlsx::read.xlsx(1) |>
  mutate(timestamp = as.POSIXct(timestamp, tz = "UTC")) |>
  select(individual_id, timestamp, x, y, buffer) |>
  as.data.table()

#### Check recaps
recaps
range(recaps$timestamp)

#### Write to file
qs::qsave(recaps, here_data("input", "mefs", "recaps.qs"))


#### End of code. 
###########################
###########################