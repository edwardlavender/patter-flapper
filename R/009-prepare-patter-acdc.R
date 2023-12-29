###########################
###########################
#### prepare-patter-acdc.R

#### Aims
# 1) Prepare ACDCPF algorithm components
#    * ACDC containers

#### Prerequisites
# 1) Obtain raw data
# 2) https://github.com/edwardlavender/flapper_appl/blob/master/R/define_global_param.R

#### Details
# In pf_forward(), for ACPF, we eliminate particles incompatible with container
# dynamics via acs_filter_container(). For the ACDCPF algorithm, we also need to
# account for depth observations to mitigate convergence issues (when an individual
# can step into the next acoustic container, but not reach valid depth cells within
# that container). For ACDCPF, we eliminate particles that are not within the required 
# distance ( mobility times the number of time steps) of at least one valid location 
# within the acoustic container. This code identifies for each unique container/depth 
# observation the set of valid cells. We will read these layers into the ACDCPF 
# algorithm to eliminate invalid proposals. 


###########################
###########################
#### Set up 

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(dv)
library(patter)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(tictoc)

#### Load data
dv::src()
bathy     <- terra::rast(here_data("spatial", "bathy.tif"))
bset      <- terra::rast(here_data("spatial", "bset.tif"))
moorings  <- readRDS(here_data("mefs", "moorings.rds"))
acoustics <- readRDS(here_data("mefs", "acoustics.rds"))
archival  <- readRDS(here_data("mefs", "archival.rds"))
ewin      <- readRasterLs(here_data("input", "depth-window"), index = FALSE)

#### Define pars
run <- TRUE


###########################
###########################
#### Identify valid cells @ detections

#### Collate data for all multiple individuals
dlist <- pat_setup_data(.acoustics = acoustics, 
                        .archival = archival, 
                        .moorings = moorings, 
                        .bathy = bathy, 
                        .lonlat = FALSE)
dlist$spatial$bset      <- bset
dlist$algorithm$ewindow <- ewin

#### Identify detection container receiver combinations 
# Build containers
containers_rc <- acs_setup_containers_rc(.dlist = dlist, .split = "individual_id")
head(containers_rc)
# Build directories for each container
if (run) {
  unlink(here_data("input", "containers"), recursive = TRUE)
  dir.create(here_data("input", "containers"))
}
folders <- here_data("input", "containers", containers_rc$receiver_key)
sapply(folders, \(folder) dir.create(folder, showWarnings = FALSE))
stopifnot(all(containers_rc$receiver_key %in% list.files(here_data("input", "containers"))))

#### Build a list of detection containers
containers <- acs_setup_detection_containers(.dlist = dlist, .rc = containers_rc)
containers[1:2]
stopifnot(all(containers_rc$receiver_key %in% names(containers)) & 
            all(names(containers) %in% containers_rc$receiver_key))
# Check which container(s) contain cells on the Digimap grid:
lapply(seq_len(length(containers)), function(i) {
  # print(i)
  container <- containers[[i]]
  if (FALSE) {
    terra::plot(container)
    terra::plot(bset, add = TRUE)
    terra::plot(container, add = TRUE)
  }
  out <- exactextractr::exact_extract(bset, sf::st_as_sf(container))[[1]]
  out |> 
    mutate(value = factor(as.character(value), levels = c("0", "1"), labels = c("bathy", "digi"))) |>
    count(value, .drop = FALSE) |> 
    tidyr::pivot_wider(names_from = value, values_from = n) |>
    mutate(index = i, container = names(containers)[i]) |> 
    select(index, container, bathy, digi) |> 
    as.data.table()
}) |> rbindlist()

#### Identify detection container receiver/depth combinations 
# Define combinations (~30 s)
tic()
containers_rcd <- acs_setup_containers_rcd(.dlist = dlist)
toc()
# Examine data
head(containers_rcd)
tail(containers_rcd)
# Generic checks 
stopifnot(length(unique(containers_rcd$index)) == 
            length(unique(containers_rcd$receiver_id_next_key))
            )
stopifnot(all(containers_rcd$receiver_id_next_key %in% containers_rc$receiver_key))
# Note that containers_rcd does not contain all receiver_keys in containers_rc
# (containers_rcd focuses on sections of overlapping time series only)
containers_rc$receiver_key[!(containers_rc$receiver_key %in% containers_rcd$receiver_id_next_key)]
# Spot checks
# * For individual 24 @ 2016-03-21 16:48:00:
# - the NEXT detection is at receiver 18
# - the NEXT depth observation at the time of the next detection is 51.19 m
acoustics[individual_id == 24, ] # |> View() 
archival[individual_id == 24, ]  # |> View() 
stopifnot(containers_rcd[timestamp == as.POSIXct("2016-03-21 16:48:00") & 
                           receiver_id_next_key == 18]$depth == 51.19)
# For individual 20 @ 2016-05-10 23:34:00
# - the NEXT detection is at receiver 27 
# - The NEXT depth observation at the time of the next detection is 17.09 m
acoustics[individual_id == 20 & timestamp > as.POSIXct("2016-05-10 23:00:00"), ] # |> View()
archival[individual_id == 20 & timestamp > as.POSIXct("2016-05-10 23:00:00")]    # |> View()
stopifnot(containers_rcd[timestamp == as.POSIXct("2016-05-10 23:34:00") & 
                           receiver_id_next_key == 27]$depth == 17.09)

#### Identify detection container cells (~8 mins on 12 forked cores)
# We write valid cells to ./data/input/containers/{container}/{depth.parquet}
gc()
tic()
if (run) {
  # .dlist = dlist
  # .containers = containers
  # .rcd = containers_rcd
  # .cl = NULL
  acs_setup_container_cells(.dlist = dlist, 
                            .containers = containers, 
                            .rcd = containers_rcd, 
                            .cl = NULL)
}
toc()

# Check file size (~500 MB)
pf_files_size(here_data("input", "containers"), recursive = TRUE)

# In pf_forward(), we account for ACDC detection container dynamics
# via acs_filter_container_acdc(). 


#### End of code. 
###########################
###########################