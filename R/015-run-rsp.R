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
# try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
dv::src()
library(ggplot2)

#### Load data 
iteration         <- qs::qread(here_data("input", "iteration", "rsp.qs"))
moorings          <- qs::qread(here_data("input", "mefs", "moorings.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))
pars              <- qs::qread(here_data("input", "pars.qs"))


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
if (FALSE) {
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
}

#### Estimate UDs
if (FALSE) {
  # Time trial (~25 s)
  lapply_estimate_ud_dbbmm(iteration = iteration[1, ], 
                           cl = NULL, 
                           plot = FALSE)
  # Implementation (7 hr, 10 cl)
  lapply_estimate_ud_dbbmm(iteration = iteration, 
                           cl = 10L, 
                           plot = FALSE)
  # (optional) Examine selected UDs
  lapply_qplot_ud(iteration, "dbbmm", "ud.tif")
}

#### Convergence (100 %)
ud_data <- lapply(file.path(iteration$folder_ud, "dbbmm", "data.qs"), function(data.qs){
  qs::qread(data.qs)
}) |> rbindlist()
table(ud_data$success)
# Verify: TRUE
file.exists(file.path(iteration$folder_ud, "dbbmm", "ud.tif")) |> table()

#### Mapping (~20 x 3 s)
if (FALSE) {
  unique(iteration$er.ad)
  length(unique(iteration$individual_id))
  length(unique(iteration$month_id))
  mapply(c(500, 250, 750), c("best", "restricted", "flexible"), 
         FUN = function(param, label) {
           # Define png args 
           png_args <- 
             list(filename = here_fig("analysis", glue("map-rsp-{label}.png")), 
                  height = 10, width = 9, units = "in", res = 800)
           # Collect data and make figure 
           iteration |> 
             filter(er.ad == param) |>
             mutate(mapfile = file.path(folder_ud, "dbbmm", "ud.tif"), 
                    individual_id = factor(individual_id, levels = sort(unique(individual_id)))) |>
             dplyr::select(row = individual_id, column = month_id, mapfile) |>
             as.data.table() |> 
             ggplot_maps(png_args = png_args)
           
         }) |> invisible()
}

#### Compute residency (~9 s)
iteration_res <- copy(iteration)
iteration_res[, file := file.path(iteration$folder_ud, "dbbmm", "ud.tif")]
iteration_res[, file_exists := file.exists(file)]
residency <- lapply_estimate_residency_ud(files = iteration_res$file[iteration_res$file_exists])
# Write output
residency <- 
  left_join(iteration_res, residency, by = "file") |> 
  mutate(algorithm = "RSP", 
         sensitivity = factor(er.ad, levels = c(500, 250, 750), labels = c("Best", "AC(-)", "AC(+)"))) |>
  select(individual_id, month_id, unit_id, algorithm, sensitivity, zone, time) |> 
  arrange(individual_id, month_id, unit_id, algorithm, sensitivity, zone) |>
  as.data.table()
qs::qsave(residency, here_data("output", "analysis-summary", "residency-rsp.qs"))


#### End of code.
###########################
###########################