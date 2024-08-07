###########################
###########################
#### run-rsp.R

#### Aims
# 1) Run RSP algorithm 

#### Prerequisites
# 1) Previous scripts


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
dv::src()

#### Load data 
iteration <- qs::qread(here_data("input", "iteration", "rsp.qs"))
moorings  <- qs::qread(here_data("input", "mefs", "moorings.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))


###########################
###########################
#### Run algorithm 

#### Set up
nrow(iteration)
moorings[, receiver_gamma := 1500]
iteration[, file_coord := file.path(folder_coord, "coord.qs")]
datasets <- list(detections_by_unit = acoustics_by_unit, moorings = moorings)

#### (optional) Testing
test <- TRUE
if (test) {
  iteration <- iteration[1:2]
}


#### Estimate coordinates
# Time trial
lapply_estimate_coord_rsp(iteration = iteration[1, ], datasets = datasets)
# Implementation (~10 mins)
lapply_estimate_coord_rsp(iteration = iteration, datasets = datasets)
# (optional) Examine selected coords
lapply_qplot_coord(iteration, 
                  "coord.qs",
                   extract_coord <- function(coord) {
                     cbind(coord$detections[[1]]$Longitude, 
                           coord$detections[[1]]$Latitude) |> 
                       terra::vect(crs = "EPSG:4326") |> 
                       terra::project("EPSG:32629") |> 
                       terra::crds() |> 
                       as.data.frame()
                   })

#### Estimate UDs
# Time trial (~25 s)
lapply_estimate_ud_dbbmm(iteration = iteration[1, ], 
                         cl = NULL, 
                         plot = FALSE)
# Implementation (~2.5 hours, 1 cl)
lapply_estimate_ud_dbbmm(iteration = iteration, 
                         cl = 2L, 
                         plot = FALSE)
# (optional) Examine selected UDs
lapply_qplot_ud(iteration, "dbbmm", "ud.tif")


#### End of code.
###########################
###########################