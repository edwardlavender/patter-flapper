###########################
###########################
#### explore-data.R

#### Aims
# 1) Examine MEFS datasets

#### Prerequisites
# 1) Process MEFS data


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)

#### Essential packages
library(dv)
library(data.table)
library(dtplyr)
library(dplyr)
library(ggplot2)
library(patter)

#### Load data
bathy     <- terra::rast(here_data("spatial", "bathy.tif"))
moorings  <- readRDS(here_data("mefs", "moorings.rds"))
acoustics <- readRDS(here_data("mefs", "acoustics.rds"))
archival  <- readRDS(here_data("mefs", "archival.rds"))


###########################
###########################
#### Examine data

# Identify individuals with archival data
acc <- 
  acoustics |> 
  filter(individual_id %in% archival$individual_id) |> 
  as.data.table()

# Count the number of observations per individual
individual_id <- 
  acc |>
  group_by(individual_id) |> 
  count(sort = TRUE) |>
  pull(individual_id)

acc$individual_id <- factor(acc$individual_id, levels = individual_id)

# Identify associated archival data
arc <- 
  archival |> 
  mutate(individual_id = factor(individual_id, levels = levels(acc$individual_id))) |>
  filter(individual_id %in% acc$individual_id) |> 
  as.data.table()

# Examine concurrent acoustic & archival time series 
ggplot() + 
  geom_line(aes(timestamp, depth * -1), data = arc) +
  geom_point(aes(timestamp, 0), data = acc) + 
  facet_wrap(~individual_id)

# Identify a suitable test individual/month with sufficient acoustic & archival data 
acc <- 
  acc |>
  filter(individual_id == 25) |> 
  filter(timestamp >= as.POSIXct("2016-07-01")) |> 
  filter(timestamp <= as.POSIXct("2016-08-01")) |> 
  # arrange(timestamp) |> 
  as.data.table()
arc <- 
  arc |> 
  filter(individual_id == acc$individual_id[1]) |> 
  # arrange(timestamp) |> 
  as.data.table()

# Calculate gaps between sequential detections
acc$gap <- Tools4ETS::serial_difference(acc$timestamp, units = "days")
plot(acc$timestamp, acc$gap, type = "l")

# (optional) Collate acoustic & archival time series 
acc$gap <- NULL
obs <- acs_setup_obs(acc, arc, 
                     .step = "2 mins", 
                     .mobility = 500, 
                     .detection_range = 750)

# Visualise time series for selected individual
obs |> 
  ggplot() + 
  geom_line(aes(timestamp, depth * -1)) + 
  geom_point(aes(timestamp, 1), data = obs |> filter(detection == 1))


#### End of code. 
###########################
###########################