###########################
###########################
#### prepare-patter-acdc.R

#### Aims
# 1) Prepare ACDC containers

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

#### Load data
dv::src()
bathy     <- terra::rast(here_data("spatial", "bathy.tif"))
bset      <- terra::rast(here_data("spatial", "bset.tif"))
moorings  <- readRDS(here_data("mefs", "moorings.rds"))
acoustics <- readRDS(here_data("mefs", "acoustics.rds"))
archival  <- readRDS(here_data("mefs", "archival.rds"))


###########################
###########################
#### Identify valid cells @ detections

#### Collate data for all multiple individuals
dlist <- pat_setup_data(.acoustics = acoustics, 
                        .archival = archival, 
                        .moorings = moorings, 
                        .bathy = bathy, 
                        .lonlat = FALSE)
dlist$spatial$bset <- bset

#### Identify detection container receiver combinations 
# Build containers
containers_rc <- acs_setup_containers_rc(.dlist = dlist, .split = "individual_id")
head(containers_rc)
# Build directories for each container
folders <- here_data("input", "containers", containers_rc$receiver_key)
sapply(folders, \(folder) dir.create(folder, showWarnings = FALSE))

#### Build a list of detection containers
containers <- acs_setup_detection_containers(.dlist = dlist, .rc = containers_rc)
containers[1:2]
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
# Define combinations 
containers_rcd <- acs_setup_containers_rcd(.dlist = dlist)
# Add index (FIX)
containers_rcd[, index := as.integer(factor(containers_rcd$receiver_id_next_key, 
                                            unique(containers_rcd$receiver_id_next_key)))]
head(containers_rcd)
tail(containers_rcd)

#### Identify detection container cells
# We write valid cells to ./data/input/containers/{container}/{depth.parquet}
acs_setup_container_cells(.dlist = dlist, 
                          .containers = containers, 
                          .rcd = containers_rcd)



