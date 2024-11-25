###########################
###########################
#### run-coa.R

#### Aims
# 1) Run COA-KUD algorithm 

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

#### Load data 
bathy             <- terra::rast(here_data("spatial", "bathy.tif"))
iteration         <- qs::qread(here_data("input", "iteration", "coa.qs"))
moorings          <- qs::qread(here_data("input", "mefs", "moorings.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))


###########################
###########################
#### Run algorithm 

#### Set up
nrow(iteration)
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
  lapply_estimate_coord_coa(iteration = iteration[1, ], datasets = datasets)
  # Implementation (~7 s)
  lapply_estimate_coord_coa(iteration = iteration, datasets = datasets)
  # (optional) Examine selected coords
  lapply_qplot_coord(iteration, "coord.qs")
}

#### Estimate UDs
if (FALSE) {
  # Time trial
  lapply_estimate_ud_spatstat(iteration = iteration[1, ], 
                              extract_coord = NULL,
                              cl = NULL, 
                              plot = FALSE)
  # Implementation (~3.78 mins: 500 pixel, 10 cl)
  lapply_estimate_ud_spatstat(iteration = iteration, 
                              extract_coord = NULL,
                              cl = 10L, 
                              plot = FALSE)
  # (optional) Examine selected UDs
  lapply_qplot_ud(iteration, "spatstat", "h", "ud.tif")
}


#### Convergence 
# 141/144 = 98 %
ud_data <- lapply(file.path(iteration$folder_ud, "spatstat", "h", "data.qs"), function(data.qs){
  qs::qread(data.qs)
}) |> rbindlist()
ud_data[!is.na(error), ]

#### Mapping (~10 s x 3)
if (FALSE) {
  unique(iteration$delta_t)
  length(unique(iteration$individual_id))
  length(unique(iteration$month_id))
  mapply(c("2 days", "1 day", "3 days"), c("best", "restricted", "flexible"), 
         FUN = function(param, label) {
           # Define png args 
           png_args <- 
             list(filename = here_fig("analysis", glue("map-coa-{label}.png")), 
                  height = 10, width = 10, units = "in", res = 600)
           # Collect data and make figure 
           iteration |> 
             filter(delta_t == param) |>
             mutate(mapfile = file.path(folder_ud, "spatstat", "h", "ud.tif"), 
                    individual_id = factor(individual_id, levels = sort(unique(individual_id)))) |>
             dplyr::select(row = individual_id, column = month_id, mapfile) |>
             as.data.table() |> 
             ggplot_maps(png_args = png_args)
           
         }) |> invisible()
}

#### Compute residency (4 s)
# (Code modified from simulate-algorithms.R)
iteration_res <- copy(iteration)
iteration_res[, file := file.path(iteration$folder_ud, "spatstat", "h", "ud.tif")]
iteration_res[, file_exists := file.exists(file)]
residency <- lapply_estimate_residency_ud(files = iteration_res$file[iteration_res$file_exists])
# Write output
residency <- 
  left_join(iteration_res, residency, by = "file") |> 
  mutate(algorithm = "COA", 
         sensitivity = factor(delta_t, levels = c("2 days", "1 day", "3 days"), labels = c("Best", "AC(-)", "AC(+)"))) |>
  select(unit_id, algorithm, sensitivity, zone, time) |> 
  arrange(unit_id, algorithm, sensitivity, zone) |>
  as.data.table()
qs::qsave(residency, here_data("output", "analysis-summary", "residency-coa.qs"))


#### End of code.
###########################
###########################