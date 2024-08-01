###########################
###########################
#### explore-data-mefs.R

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
dv::clear()

#### Essential packages
library(dv)
library(data.table)
library(dtplyr)
library(dplyr)
library(ggplot2)
library(patter)
library(tictoc)
dv::src()

#### Load data
bathy     <- terra::rast(here_data("spatial", "bathy.tif"))
moorings  <- qs::qread(here_data("input", "mefs", "moorings.qs"))
acoustics <- qs::qread(here_data("input", "mefs", "acoustics.qs"))
archival  <- qs::qread(here_data("input", "mefs", "archival.qs"))


###########################
###########################
#### Summary statistics

# 39 individuals with data 
length(unique(c(acoustics$individual_id, archival$individual_id)))

# 33 acoustic time series 
length(unique(acoustics$individual_id))

# 21 archival time series
length(unique(archival$individual_id))

# 15/33 acoustic time series are associated with depth time series
table(unique(acoustics$individual_id) %in% unique(archival$individual_id))

# 15/21 archival time series are associated with acoustic data
table(unique(archival$individual_id) %in% unique(acoustics$individual_id))


###########################
###########################
#### Visualise time series

# ~24 s
tic()
png(here_fig("mefs", "all-time-series.png"), 
    height = 10, width = 20, res = 600, units = "in")
ggplot() +
  geom_line(aes(timestamp, depth*-1), lwd = 0.25, data = archival) + 
  geom_point(aes(timestamp, 0), colour = "red", data = acoustics) + 
  scale_x_datetime(labels = scales::date_format("%m"), 
                   breaks = "1 month") +
  facet_wrap(~individual_id, scales = "free") |> 
  print()
dev.off()
toc()


###########################
###########################
#### Define example individual

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
acc$gap <- serial_difference(acc$timestamp, units = "days")
plot(acc$timestamp, acc$gap, type = "l")

# (optional) Collate acoustic & archival time series 
# & visualise time series for a selected individual
# (TO DO)


#### End of code. 
###########################
###########################