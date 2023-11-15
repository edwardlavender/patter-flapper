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


###########################
###########################
#### Save datasets

saveRDS(acoustics, here_data("mefs", "acoustics.rds"))
saveRDS(archival, here_data("mefs", "archival.rds"))
saveRDS(moorings, here_data("mefs", "moorings.rds"))


#### End of code. 
###########################
###########################