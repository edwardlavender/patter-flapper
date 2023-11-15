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

#### Essential packages
library(dv)
library(MEFS)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(sf)

#### Load data
acoustics <- as.data.table(acoustics)
archival  <- as.data.table(archival)
moorings  <- as.data.table(moorings)
bathy <- terra::rast(here_data("spatial", "bathy.tif"))


###########################
###########################
#### Moorings

#### Enforce {patter} requirements
moorings <- 
  moorings |> 
  select(receiver_id, 
         receiver_lon = long_receiver, receiver_lat = lat_receiver, 
         receiver_start = date_operation_start_receiver, receiver_end = date_operation_end_receiver) |> 
  mutate(receiver_range = 750) |> 
  as.data.table()

#### Define UTM coordinates
rxy <- 
  moorings |>
  st_as_sf(coords = c("receiver_lon", "receiver_lat")) |> 
  st_set_crs(4326) |> 
  st_transform(terra::crs(bathy)) |> 
  st_coordinates()
moorings$receiver_easting  <- rxy[, 1]
moorings$receiver_northing <- rxy[, 2]

#### Check validity on `bathy`
# * On Howe et al. data, 7 receivers are invalid on `bathy` e.g., along lismore
# * This has been resolved by merging Howe et al & Digimap data
moorings$bathy <- terra::extract(bathy, moorings[, .(receiver_easting, receiver_northing)])$depth
table(is.na(moorings$bathy))
if (FALSE) {
  # Plot receiver positions
  terra::plot(bathy)
  text(moorings$receiver_easting, moorings$receiver_northing, moorings$receiver_id, cex = 0.5)
  # Examine bathymetry around selected receiver 
  xy <- cbind(moorings$receiver_easting[1], moorings$receiver_northing[1])
  bathy |> 
    terra::crop(terra::ext(xy[1] - 1e3, xy[1] + 1e3, xy[2] - 1e3, xy[2] + 1e3)) |> 
    terra::plot()
  points(xy)
}


###########################
###########################
#### Save datasets

saveRDS(acoustics, here_data("mefs", "acoustics.rds"))
saveRDS(archival, here_data("mefs", "archival.rds"))
saveRDS(moorings, here_data("mefs", "moorings.rds"))


#### End of code. 
###########################
###########################