###########################
###########################
#### process-data-mefs.R

#### Aims
# 1) Process MEFS data for analysis 

#### Prerequisites
# 1) Obtain raw data


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(MEFS)
library(sf)
dv::src()

#### Load data
skateids  <- as.data.table(skateids)
acoustics <- as.data.table(acoustics)
archival  <- as.data.table(archival)
moorings  <- as.data.table(moorings)
bathy     <- terra::rast(here_data("spatial", "bathy.tif"))
coast     <- qreadvect(here_data("spatial", "coast.qs"))


###########################
###########################
#### Electonic tagging & tracking datasets

#### Enforce {patter} requirements
moorings <- 
  moorings |> 
  select(receiver_id, 
         receiver_lon = long_receiver, 
         receiver_lat = lat_receiver, 
         receiver_start = date_operation_start_receiver, 
         receiver_end = date_operation_end_receiver) |> 
  mutate(receiver_start = date_to_POSIXct(receiver_start), 
         receiver_end = date_to_POSIXct(receiver_end)) |>
  as.data.table()

#### Define UTM coordinates
rxy <- 
  moorings |>
  st_as_sf(coords = c("receiver_lon", "receiver_lat")) |> 
  st_set_crs(4326) |> 
  st_transform(terra::crs(bathy)) |> 
  st_coordinates()
moorings <- 
  moorings |> 
  select(receiver_id, receiver_start, receiver_end) |> 
  mutate(receiver_x = rxy[, 1], 
         receiver_y = rxy[, 2]) |> 
  as.data.table()

#### Check validity on `bathy`
moorings$map_value <- 
  terra::extract(bathy, cbind(moorings$receiver_x, moorings$receiver_y))[, 1]
# Receiver 2 is in very shallow water (NA on high-resolution bathymetry raster)
moorings[is.na(map_value), ]
# Further examination:
if (FALSE) {
  # Plot receiver positions
  terra::plot(bathy)
  text(moorings$receiver_x, moorings$receiver_y,
       moorings$receiver_id, cex = 0.5)
  
  # Examine bathymetry around selected receiver 
  xy <- cbind(moorings$receiver_x[1], moorings$receiver_y[1])
  buf <- 1e2
  bathy |> 
    terra::crop(terra::ext(xy[1] - buf, xy[1] + buf, 
                           xy[2] - buf, xy[2] + buf)) |> 
    terra::plot()
  points(xy)
  terra::lines(coast)
  # Check depth in nearest location (use buf = 1e2 and locator())
  terra::extract(bathy, cbind(719462, 6273481))

  # There were never any detections at this receiver
  2 %in% MEFS::acoustics$receiver_id 
}
# Exclude receiver 2
# > It is incompatible with our bathymetry layer
# > (This is no longer implemented as it is compatible with the aggregated raster)
# moorings <- moorings[receiver_id != 2, ]

#### Validate datasets
dlist <- pat_setup_data(.map = bathy, 
                        .detections = acoustics, 
                        .moorings = moorings, 
                        .archival = archival)


###########################
###########################
#### Save datasets

qs::qsave(skateids, here_data("input", "mefs", "skateids.qs"))
qs::qsave(dlist$acoustics, here_data("input", "mefs", "acoustics.qs"))
qs::qsave(dlist$archival, here_data("input", "mefs", "archival.qs"))
qs::qsave(dlist$moorings, here_data("input", "mefs", "moorings.qs"))


#### End of code. 
###########################
###########################