###########################
###########################
#### prepare-patter-obs.R

#### Aims
# 1) Prepare acoustic (and associated archival) observations for analysis

#### Prerequisites
# 1) Obtain raw data

#### TO DO
# * Publication-quality plots
# * Create an abacus plot of the raw acoustic time series for all individuals
# * Include the depth time series (?) 
# * Add blocks delineating which months were selected for analysis


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
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(patter)
library(tictoc)

#### Load data 
source(here_r("002-define-helpers.R"))
skateids  <- readRDS(here_data("mefs", "skateids.rds"))
acoustics <- readRDS(here_data("mefs", "acoustics.rds"))
archival  <- readRDS(here_data("mefs", "archival.rds"))
bathy     <- terra::rast(here_data("spatial", "bathy.tif"))


###########################
###########################
#### Global processing 

#### Define study period 
start <- as.Date("2016-04-01")
end   <- as.Date("2017-05-31")
mons  <- seq(start, end, by = "months")

#### Globally process skateids
skateids$date_removal_tag_acc[is.na(skateids$date_removal_tag_acc)] <- end
skateids <- 
  skateids |> 
  mutate(arc_start_date = date_deployment_tag) |>
  rename(acc_start_date = date_deployment_tag, 
         acc_end_date = date_removal_tag_acc, 
         arc_end_date = date_removal_tag_dst
         )

#### Globally process acoustics data 
acoustics <- 
  acoustics |> 
  lazy_dt() |>
  # Restrict observations within time frame 
  filter(timestamp >= start & timestamp <= end) |>
  # Define time blocks 
  mutate(block = lubridate::floor_date(timestamp, "week"), 
         mmyy = Tools4ETS::mmyy(timestamp)) |> 
  # For each individual, focus on months entirely at liberty 
  left_join(skateids |> 
              select(individual_id, 
                     acc_start_date, 
                     acc_end_date) |> 
              mutate(acc_start_date = 
                       lubridate::ceiling_date(acc_start_date, "month"),
                     acc_end_date = 
                       lubridate::floor_date(acc_end_date, "month")) |> 
              as.data.table(),
            by = "individual_id") |>
  group_by(individual_id) |> 
  filter(timestamp >= acc_start_date & timestamp <= acc_end_date) |>
  ungroup() |>
  as.data.table()

#### (optional) Exclude individuals that moved beyond the study area
# * We exclude individuals known to have moved beyond study area
# * ... based on depth time series
# * Other individuals may have moved beyond the study area 
# * We validate algorithm convergence using those individuals
if (FALSE) {
  # The study area spans the full range of observed depths (~15 s)
  tic()
  max(archival$depth)
  terra::plot(bathy > max(archival$depth))
  terra::global(bathy, "max", na.rm = TRUE)
  toc()
}

#### Globally process archival data
archival <-
  archival |> 
  mutate(mmyy = Tools4ETS::mmyy(timestamp)) |> 
  as.data.table()


###########################
###########################
#### Identify individuals/months for analysis 

# For each individual/month, count the number of blocks (e.g. weeks) with detections
smr <- 
  acoustics |> 
  group_by(individual_id, mmyy) |>
  summarise(blocks = n_distinct(block)) |> 
  ungroup() |>
  as.data.table()

# Identify individuals & months with 'sufficient' observations
# E.g., observation(s) in at least one week, two weeks etc.
crit <- 2L
smr |> 
  filter(blocks >= crit) |>
  as.data.table()

# Count the number of individuals per mmyy category with sufficient observations
smr |> 
  filter(blocks >= crit) |>
  group_by(mmyy) |>
  arrange(individual_id, .by_group = TRUE) |>
  summarise(n = n(), 
            id_list = list(individual_id), 
            id_str = paste(individual_id, collapse = ", ")) |> 
  as.data.table()

# Define selection criterial for acoustic time series
sel <- 
  smr |> 
  filter(blocks >= crit) |>
  as.data.table()

# Collect acoustic and archival data
data_ls <- 
  lapply(split(sel, collapse::seq_row(sel)), function(d) {
    
    # Define acoustic data 
    acc <- 
      acoustics |> 
      filter(individual_id == d$individual_id[1]) |> 
      filter(mmyy %in% d$mmyy) |> 
      as.data.table()
    
    # Define archival data 
    period <- 
    arc <- 
      archival |> 
      filter(individual_id == d$individual_id[1]) |> 
      filter(mmyy %in% d$mmyy) |> 
      as.data.table()
    
    if (nrow(arc) == 0L) {
      arc <- NULL
    } else {
      # (optional) Only retain archival data that spans the whole period 
      # > Trial with/without to compare number & quality of time series for analysis
      period <- mmyyrng(d$mmyy)
      if (min(arc$timestamp) > min(period)) {
        arc <- NULL
      }
      if (max(arc$timestamp) < max(period)) {
        arc <- NULL
      }
    }
    
    # Collate data
    list(acoustics = acc, archival = arc)
  }) 



###########################
###########################
#### Collate observations

obs_ls <- 
  lapply(data_ls, function(d) {
    
    # d <- data_ls[[1]]
    
    # Define period
    period <- mmyyrng(d$acoustics$mmyy[1])
    
    # ACPF observations
    # TO DO
    # * Fix start time & end time to start/beginning of month
    acpf <- acs_setup_obs(.acoustics = d$acoustics, 
                          .step = "2 mins",
                          .period = period,
                          .mobility = 500, 
                          .detection_range = 750)
    acpf[, individual_id := d$acoustics$individual_id[1]]
    acpf[, block := d$acoustics$mmyy[1]]
    acpf[, algorithm := "acpf"]

    # ACDCPF observations
    acdcpf <- NULL
    if (!is.null(d$archival)) {
      acdcpf <- acs_setup_obs(.acoustics = d$acoustics, 
                              .archival = d$archival, 
                              .trim = FALSE,
                              .step = "2 mins", 
                              .period = period,
                              .mobility = 500, 
                              .detection_range = 750)
      acdcpf[, individual_id := d$acoustics$individual_id[1]]
      acdcpf[, block := d$acoustics$mmyy[1]]
      acdcpf[, algorithm := "acdcpf"]
    }
    
    # Collate observations
    list(acpf = acpf, 
         acdcpf = acdcpf)
  })

#### Visualise observations
# ACPF
tic()
png(here_fig("all-obs-acpf.png"),
    height = 10, width = 12, units = "in", res = 600)
obs_ls |> 
  purrr::list_flatten() |> 
  rbindlist(fill = TRUE) |> 
  filter(algorithm == "acpf") |>
  ggplot() + 
  geom_point(aes(timestamp, factor(individual_id)), data = . %>% filter(detection == 1L)) + 
  facet_wrap(~block, scales = "free")
dev.off()
# ACDCPF
png(here_fig("all-obs-acdcpf.png"),
    height = 10, width = 12, units = "in", res = 600)
obs_ls |> 
  purrr::list_flatten() |> 
  rbindlist(fill = TRUE) |> 
  filter(algorithm == "acdcpf") |>
  ggplot() + 
  geom_line(aes(timestamp, depth * -1), lwd = 0.5) +
  geom_point(aes(timestamp, 0), 
             data = . %>% filter(detection == 1L), 
             colour = "red") + 
  scale_x_datetime(labels = scales::date_format("%d")) +
  facet_wrap(~block + individual_id, scales = "free")
dev.off()
toc()


#### End of code. 
###########################
###########################