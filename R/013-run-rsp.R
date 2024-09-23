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
library(ggplot2)
dv::src()

#### Load data 
iteration <- qs::qread(here_data("input", "iteration", "rsp.qs"))
moorings  <- qs::qread(here_data("input", "mefs", "moorings.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))
pars      <- qs::qread(here_data("input", "pars.qs"))


###########################
###########################
#### Run algorithm 

#### Set up
nrow(iteration)
moorings[, receiver_gamma := pars$pdetection$receiver_gamma[1]]
iteration[, file_coord := file.path(folder_coord, "coord.qs")]
datasets <- list(detections_by_unit = acoustics_by_unit, moorings = moorings)

#### (optional) Testing
test <- FALSE
if (test) {
  iteration <- iteration[1:2]
}

#### Estimate coordinates
# Time trial
lapply_estimate_coord_rsp(iteration = iteration[1, ], datasets = datasets)
# Implementation (~6 mins)
lapply_estimate_coord_rsp(iteration = iteration, datasets = datasets)
# (optional) Examine selected coords
lapply_qplot_coord(iteration, 
                  "coord.qs",
                   extract_coord = function(coord) {
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
# Implementation (26 mins, 8 cl; 15 mins, 11 cl)
lapply_estimate_ud_dbbmm(iteration = iteration, 
                         cl = 11L, 
                         plot = FALSE)
# (optional) Examine selected UDs
lapply_qplot_ud(iteration, "dbbmm", "ud.tif")

#### Mapping (~10 x 3 s)
unique(iteration$er.ad)
length(unique(iteration$individual_id))
length(unique(iteration$month_id))
mapply(c(100, 50, 150), c("best", "restricted", "flexible"), 
       FUN = function(param, label) {
         # Define png args 
         png_args <- 
           list(filename = here_fig("analysis", glue("map-rsp-{label}.png")), 
                height = 5, width = 10, units = "in", res = 600)
         # Collect data and make figure 
         iteration |> 
           filter(er.ad == param) |>
           mutate(mapfile = file.path(folder_ud, "dbbmm", "ud.tif"), 
                  individual_id = factor(individual_id, levels = sort(unique(individual_id)))) |>
           select(row = individual_id, column = month_id, mapfile) |>
           as.data.table() |> 
           ggplot_maps(png_args = png_args)
         
       }) |> invisible()


#### End of code.
###########################
###########################