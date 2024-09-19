###########################
###########################
#### process-data-mefs.R

#### Aims
# 1) Prepare acoustic (and associated archival) observations for analysis

#### Prerequisites
# 1) Obtain raw data


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
# try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(ggplot2)
dv::src()

#### Load data 
skateids  <- qs::qread(here_data("input", "mefs", "skateids.qs"))
moorings  <- qs::qread(here_data("input", "mefs", "moorings.qs"))
acoustics <- qs::qread(here_data("input", "mefs", "acoustics.qs"))
archival  <- qs::qread(here_data("input", "mefs", "archival.qs"))
bathy     <- terra::rast(here_data("spatial", "bathy.tif"))


###########################
###########################
#### Global data processing

###########################
#### Skate IDs

start       <- as.Date("2016-04-01")
end         <- as.Date("2017-05-31")
mons        <- seq(start, end, by = "months")

skateids$date_removal_tag_acc[is.na(skateids$date_removal_tag_acc)] <- end
skateids <- 
  skateids |> 
  mutate(arc_start_date = date_deployment_tag) |>
  rename(acc_start_date = date_deployment_tag, 
         acc_end_date = date_removal_tag_acc, 
         arc_end_date = date_removal_tag_dst
  )


###########################
#### Acoustics

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
#### Archival

#### Inclusion criteria
# In each month, we include individuals that:
# * were at liberty for the entire month
# * have complete depth time series for that month 
# (for the reasons described above)

#### Processing
# Round time stamps
archival[, timestamp := lubridate:::round_date(timestamp, "2 mins")]
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
# Additional check for NAs
any(is.na(archival$depth))

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
#### Define iteration dataset

#### (optional) Focus on individuals/months with both acoustic & archival data
# We do not analyse all data
# We focus on individuals/months with both datasets
# This facilitates comparisons between DCPF, ACPF and ACDCPF algorithms 
acoustics[, key := paste(individual_id, mmyy)]
archival[, key := paste(individual_id, mmyy)]
keys      <- intersect(acoustics$key, archival$key)
acoustics <- acoustics[key %in% keys, ]
archival  <- archival[key %in% keys, ]
# There are 50 individual/month combinations
length(keys)
# For each individual, this is the number of months with data:
acoustics |> group_by(individual_id) |> summarise(n_months = length(unique(mmyy)))
identical(acoustics |> group_by(individual_id) |> summarise(n_months = length(unique(mmyy))),
          archival |> group_by(individual_id) |> summarise(n_months = length(unique(mmyy)))) |> stopifnot()
# For each month, this is the number of individuals with data
acoustics |> group_by(mmyy) |> summarise(n_ids = length(unique(individual_id)))
identical(acoustics |> group_by(mmyy) |> summarise(n_ids = length(unique(individual_id))) ,
          archival |> group_by(mmyy) |> summarise(n_ids = length(unique(individual_id)))) |> stopifnot()
# Update mons vector (not required)
# mons <- mons[mons %in% as.Date(archival$timestamp)]

#### Define template individual/months dataset
# Individuals/months are the unit of analysis
units <- 
  CJ(individual_id = sort(unique(c(acoustics$individual_id, archival$individual_id))),
     mmyy = mmyy(mons)) |> 
  arrange(individual_id) |> 
  as.data.table()
units[, unit_id := seq_len(.N)]

#### Build a list of raw acoustic datasets
acoustics_by_unit <- 
  lapply(split(units, units$unit_id), function(d) {
  
  # Identify acoustic data
  # d = units[1, ]
  acc <- 
    acoustics |> 
    filter(individual_id == d$individual_id[1]) |> 
    filter(mmyy %in% d$mmyy) |> 
    as.data.table()
  
  if (nrow(acc) == 0L) {
    return(NULL)
  }
  
  # Retain data for individuals at liberty for the entire month
  # > Achieved above.
  
  # Retain data for individuals with detections in > two weeks on a total of >= 7 days
  # > Achieved above.
  
  acc

})

#### Build a list of raw archival datasets
archival_by_unit <- 
  lapply(split(units, units$unit_id), function(d) {
  
  # Identify archival data
    arc <- 
      archival |> 
      filter(individual_id == d$individual_id[1]) |> 
      filter(mmyy %in% d$mmyy) |> 
      as.data.table()
    
    if (nrow(arc) == 0L) {
      return(NULL)
    }
    
  # Retain data for individuals at liberty for the entire month
  # > Achieved above.
  
  # Retain data for individuals with complete archival time series
  # > Achieved above. 

  arc

})

#### Link units (individuals/months) to data availability
unitsets <- lapply(split(units, units$unit_id), function(d) {
  
  acc <- acoustics_by_unit[[d$unit_id]]
  arc <- archival_by_unit[[d$unit_id]]
  
  if (is.null(acc) & is.null(arc)) {
    return(NULL)
  }
  if (!is.null(acc) & is.null(arc)) {
    datasets <- "ac"
  }
  if (is.null(acc) & !is.null(arc)) {
    datasets <- "dc"
  }
  if (!is.null(acc) & !is.null(arc)) {
    datasets <- c("ac", "dc", "acdc")
  }
  
  lapply(datasets, function(dset) {
    data.table(unit_id = d$unit_id,
               individual_id = d$individual_id, 
               month_id = d$mmyy,
               dataset = dset)
  }) |> rbindlist()
  
}) |> rbindlist()

head(unitsets)
table(unitsets$dataset)

#### Checks
# For every individual/month, there is an ac, dc and acdc dataset
unitsets |> 
  group_by(individual_id, month_id) |> 
  summarise(n = n()) |> 
  pull(n) |> 
  table()


###########################
###########################
#### Visualise time series

# Visualise ACPF time series

# Visualise DCPF time series

# Visualise ACDCPF time series 


###########################
###########################
#### Save outputs

# ~3 s
tic()
qs::qsave(acoustics_by_unit, here_data("input", "acoustics_by_unit.qs"))
qs::qsave(archival_by_unit, here_data("input", "archival_by_unit.qs"))
qs::qsave(unitsets, here_data("input", "unitsets.qs"))
toc()


#### End of code. 
###########################
###########################