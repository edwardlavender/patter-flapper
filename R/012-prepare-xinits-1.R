###########################
###########################
#### prepare-xinits-1.R

#### Aims
# (1) Prepare initial locations for the filter (part 1)
# * We use recapture events & acoustic detections to define xinit
# * This code prepares the recaptures dataset

#### Prerequisites
# (1) Download of skatespotter dataset (Lacie_Share)


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
map       <- terra::rast(here_data("spatial", "bathy.tif"))
iteration <- qs::qread(here_data("input", "iteration", "patter.qs"))
skateids  <- qs::qread(here_data("input", "mefs", "skateids.qs"))
recaps    <- readRDS(here_data_raw("movement", "recaptures_processed.rds"))


###########################
###########################
#### Prepare recaps

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
    recaps_skatespotter |> select("individual_id", "timestamp", "x", "y"),
    ignore.attr = TRUE
  )

# Save
qs::qsave(recaps, here_data("input", "mefs", "recaps.qs"))


#### End of code. 
###########################
###########################