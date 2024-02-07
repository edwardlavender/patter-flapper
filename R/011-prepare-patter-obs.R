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
dv::src()
pars      <- readRDS(here_data("input", "pars.rds"))
skateids  <- readRDS(here_data("mefs", "skateids.rds"))
acoustics <- readRDS(here_data("mefs", "acoustics.rds"))
archival  <- readRDS(here_data("mefs", "archival.rds"))
bathy     <- terra::rast(here_data("spatial", "bathy.tif"))


###########################
###########################
#### Define parameters

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


###########################
###########################
#### Process acoustics

#### Inclusion criteria
# In each month, we include individuals that:
# * were at liberty for the entire month
# * were detected in at least two different 7 day periods
# The liberty criterion ensures comparability among individual maps
# The detection criterion excludes individuals with the poorest time series
# (e.g., individuals that were detected following tagging and then appeared to leave the array)

#### Initial processing  
acoustics <- 
  acoustics |> 
  lazy_dt() |>
  # Restrict observations within time frame 
  mutate(date = as.Date(timestamp)) |>
  filter(date >= start & date <= end) |>
  # Define time blocks 
  mutate(block = lubridate::floor_date(timestamp, "week"), 
         mmyy = mmyy(timestamp)) |>
  # For each individual, focus on months entirely at liberty 
  left_join(skateids |> 
              select(individual_id, 
                     acc_start_date, 
                     acc_end_date) |> 
              mutate(acc_start_date = 
                       as.Date(lubridate::ceiling_date(acc_start_date, "month")),
                     acc_end_date = 
                       as.Date((lubridate::floor_date(acc_end_date, "month") + lubridate::days(1)) - lubridate::days(2))) |> 
              as.data.table(),
            by = "individual_id") |>
  group_by(individual_id) |> 
  filter(date >= acc_start_date & date <= acc_end_date) |>
  ungroup() |>
  as.data.table()

#### Focus on individuals with a minimum number of detections
# For each individual/month, count the number of blocks (e.g. weeks) with detections
smr <- 
  acoustics |> 
  group_by(individual_id, mmyy) |>
  summarise(blocks = n_distinct(block), 
            days = difftime(max(timestamp), min(timestamp), units = "days")) |> 
  ungroup() |>
  as.data.table()
# Identify individuals & months with 'sufficient' observations
# E.g., observation(s) in at least one week, two weeks etc.
pass <- 
  smr |> 
  filter(blocks >= 2L & days > 7L) |>
  mutate(indicator = 1L) |>
  as.data.table()
# Count the number of individuals per mmyy category with sufficient observations
pass |>
  group_by(mmyy) |>
  arrange(individual_id, .by_group = TRUE) |>
  summarise(n = n(), 
            id_list = list(individual_id), 
            id_str = paste(individual_id, collapse = ", ")) |> 
  as.data.table()
# Filter acoustics accordingly 
nrow(acoustics)
acoustics <- 
  acoustics |> 
  left_join(pass |> 
              select(individual_id, mmyy, indicator), 
            by = c("individual_id", "mmyy")) |> 
  filter(!is.na(indicator)) |> 
  select(individual_id, mmyy, timestamp, receiver_id) |>
  arrange(individual_id, mmyy, timestamp, receiver_id) |>
  as.data.table()
nrow(acoustics)


###########################
###########################
#### Process archival data

#### Inclusion criteria
# In each month, we include individuals that:
# * were at liberty for the entire month
# * have complete depth time series for that month 
# (for the reasons described above)

#### Processing
# Round time stamps
archival[, timestamp := lubridate:::round_date(timestamp, pars$patter$step)]
stopifnot(all(archival |>
                lazy_dt() |>
                group_by(individual_id, timestamp) |> 
                summarise(n = n()) |> 
                pull(n) == 1L))
# Continue processing
archival <- 
  archival |> 
  lazy_dt() |>
  # Restrict observations within time frame 
  mutate(date = as.Date(timestamp)) |>
  filter(date >= start & date <= end) |>
  # Define time blocks 
  mutate(block = lubridate::floor_date(timestamp, "week"), 
         mmyy = mmyy(timestamp)) |> 
  # For each individual, focus on months entirely at liberty 
  left_join(skateids |> 
              select(individual_id, 
                     arc_start_date, 
                     arc_end_date) |> 
              mutate(arc_start_date = 
                       as.Date(lubridate::ceiling_date(arc_start_date, "month")),
                     arc_end_date = 
                       as.Date((lubridate::floor_date(arc_end_date, "month") + lubridate::days(1)) - lubridate::days(2))) |> 
              as.data.table(),
            by = "individual_id") |>
  group_by(individual_id) |> 
  filter(date >= arc_start_date & date <= arc_end_date) |>
  ungroup() |>
  as.data.table()

#### (optional) Exclude individuals that moved beyond the study area
# * We exclude individuals known to have moved beyond study area based on depth time series
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


###########################
###########################
#### Identify individuals/months for analysis 

combs <- 
  CJ(
    individual_id = unique(c(acoustics$individual_id, archival$individual_id)), 
    mmyy = mmyy(mons)) |> 
  arrange(individual_id)

data_ls <- 
  pbapply::pblapply(split(combs, collapse::seq_row(combs)), function(d) {
    
    # Define acoustic data 
    acc <- 
      acoustics |> 
      filter(individual_id == d$individual_id[1]) |> 
      filter(mmyy %in% d$mmyy) |> 
      as.data.table()
    if (nrow(acc) == 0L) {
      acc <- NULL
    }
    
    # Define archival data 
    arc <- 
      archival |> 
      filter(individual_id == d$individual_id[1]) |> 
      filter(mmyy %in% d$mmyy) |> 
      as.data.table()
    if (nrow(arc) == 0L) {
      arc <- NULL
    }
    
    # Collate data
    if (is.null(acc) && is.null(arc)) {
      return(NULL)
    }
    list(acoustics = copy(acc), archival = copy(arc))
    
  }) |> 
  plyr::compact()


###########################
###########################
#### Collate observations

#### Collate observations 
obs_ls <- 
  lapply(seq_len(length(data_ls)), function(i) {
    
    print(i)
    d <- data_ls[[i]]
    
    # Define period
    if (!is.null(d$acoustics)) {
      period <- mmyyrng(d$acoustics$mmyy[1])
    } else if (!is.null(d$archival)) {
      period <- mmyyrng(d$archival$mmyy[1])
    } 

    # DCPF observations
    dcpf <- NULL
    if (!is.null(d$archival)) {
      dcpf <- acs_setup_obs(.archival = d$archival, 
                            .step = pars$patter$step,
                            .period = period, 
                            .mobility = pars$patter$mobility)
      # Exclude individuals with gaps (due to recapture events)
      if (any(is.na(dcpf$depth))) {
        dcpf <- NULL
      } else {
        dcpf[, individual_id := d$archival$individual_id[1]]
        dcpf[, block := as.character(d$archival$mmyy[1])]
        dcpf[, algorithm := "dcpf"]
      }
    }
    
    # ACPF observations
    acpf <- NULL
    if (!is.null(d$acoustics)) {
      acpf <- acs_setup_obs(.acoustics = d$acoustics, 
                            .step = pars$patter$step,
                            .period = period,
                            .mobility = pars$patter$mobility, 
                            .detection_range = pars$patter$detection_range)
      acpf[, individual_id := d$acoustics$individual_id[1]]
      acpf[, block := as.character(d$acoustics$mmyy[1])]
      acpf[, algorithm := "acpf"]
    }

    # ACDCPF observations
    acdcpf <- NULL
    if (!is.null(d$acoustics) && !is.null(d$archival)) {
      acdcpf <- acs_setup_obs(.acoustics = d$acoustics, 
                              .archival = d$archival, 
                              .trim = FALSE,
                              .step = pars$patter$step, 
                              .period = period,
                              .mobility = pars$patter$mobility, 
                              .detection_range = pars$patter$detection_range)
      acdcpf[, individual_id := d$acoustics$individual_id[1]]
      acdcpf[, block := as.character(d$acoustics$mmyy[1])]
      acdcpf[, algorithm := "acdcpf"]
    }
    
    # Collate observations
    list(acpf = acpf, 
         dcpf = dcpf, 
         acdcpf = acdcpf)
  }) 
# Define dataframe
obs_data <-   
  obs_ls |> 
  purrr::list_flatten() |> 
  rbindlist(fill = TRUE)

#### Validation
# Confirm that the first/last elements are the first/last day in each month 
val <- 
  obs_data |> 
  group_by(individual_id, block, algorithm) |> 
  summarise(first = timestamp[1], 
            last = timestamp[n()]) |> 
  ungroup()
stopifnot(all(as.Date(val$first) == 
                lubridate::floor_date(val$first, "month")))
stopifnot(all(as.Date(val$last) == 
                lubridate::ceiling_date(val$last, "month") - lubridate::days(1)))


###########################
###########################
#### Include directories 

#### Define directories
# data/output/forward/{individual}/{block}/{algorithm}/output/
# > log.txt
# > history/
# > diagnostics 

#### Define folders
folders <- 
  obs_data |>
  select(individual_id, block, algorithm) |>
  group_by(individual_id, block, algorithm) |>
  slice(1L) |>
  mutate(folder = here_data("output", "forward", individual_id, block, algorithm, "output")) |>
  as.data.table()

#### Include in obs_data
obs_data <- 
  obs_data |>
  left_join(folders, by = c("individual_id", "block", "algorithm")) |>
  as.data.table()
head(obs_data)

#### Define single-level list for algorithm implementations
obs_ls <- split(obs_data, by = c("individual_id", "block", "algorithm"), drop = TRUE)
stopifnot(!any(sapply(obs_ls, nrow) == 0L))
str(obs_ls[[1]])
length(obs_ls)


###########################
###########################
#### Visualise observations

# ~38 s

if (FALSE) {
  
  # ACPF
  tic()
  png(here_fig("all-obs-acpf.png"),
      height = 10, width = 12, units = "in", res = 600)
  obs_data |> 
    filter(algorithm == "acpf") |>
    ggplot() + 
    geom_point(aes(timestamp, factor(individual_id)), data = . %>% filter(detection == 1L)) + 
    facet_wrap(~block, scales = "free")
  dev.off()
  
  # DCPF
  png(here_fig("all-obs-dcpf.png"),
      height = 10, width = 12, units = "in", res = 600)
  obs_data |> 
    filter(algorithm == "dcpf") |>
    ggplot() + 
    geom_line(aes(timestamp, depth * -1), lwd = 0.5) +
    facet_wrap(~individual_id + block, scales = "free")
  dev.off()
  
  # ACDCPF
  png(here_fig("all-obs-acdcpf.png"),
      height = 10, width = 12, units = "in", res = 600)
  obs_data |> 
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
  
}


###########################
###########################
#### Save outputs

# ~3 s
tic()
qs::qsave(obs_data, here_data("input", "obs_data.qs"))
qs::qsave(obs_ls, here_data("input", "obs.qs"))
toc()


#### End of code. 
###########################
###########################