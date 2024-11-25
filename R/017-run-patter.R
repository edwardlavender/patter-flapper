###########################
###########################
#### run-patter.R

#### Aims
# 1) Run patter algorithms

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
if (!patter:::os_linux()) {
  map <- terra::rast(here_data("spatial", "bathy.tif"))
}
pars              <- qs::qread(here_data("input", "pars.qs"))
iteration         <- qs::qread(here_data("input", "iteration", "patter.qs"))
moorings          <- qs::qread(here_data("input", "mefs", "moorings.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))
archival_by_unit  <- qs::qread(here_data("input", "archival_by_unit.qs"))
behaviour_by_unit <- qs::qread(here_data("input", "behaviour_by_unit.qs"))


###########################
###########################
#### Run algorithm 

#### Julia Set up
if (patter:::os_linux()) {
  stopifnot(!any(c("terra", "sf") %in% sort(loadedNamespaces())))
}
julia_connect()
set_seed()
set_map(here_data("spatial", "bathy.tif"))
set_model_move_components()

#### Clean up
if (FALSE) {
  unlink(iteration$folder_coord, recursive = TRUE)
  unlink(iteration$folder_ud, recursive = TRUE)
  list.files(iteration$folder_coord, recursive = TRUE)
  list.files(iteration$folder_ud, recursive = TRUE)
  dirs.create(iteration$folder_coord)
  dirs.create(iteration$folder_ud)
}

#### Prepare iterations/datasets
iteration[, folder_xinit := file.path("data", "input", "xinit", "1.5", individual_id, month_id)]
iteration[, file_coord := file.path(folder_coord, "coord-smo.qs")]
datasets <- list(detections_by_unit = acoustics_by_unit, 
                 moorings           = moorings,
                 archival_by_unit   = archival_by_unit, 
                 behaviour_by_unit  = behaviour_by_unit)

#### Review inputs
# Check iterations
# * Each unit_id is a individual/month combination
# * For each unit_id, we have AC/DC/ACDC algorithm implementations
# * For each algorithm implementation, we have different parameterisations (sensitivity)
head(iteration)
# Check datasets 
# * unit_ids were defined as _all possible_ individual-month combinations
# * (see process-data-mefs.R)
# -> We have 154 possible unit ids (individual-month) combinations
ids  <- unique(iteration$individual_id)
mons <- mmyy(seq(as.Date("2016-04-01"), as.Date("2017-05-31"), by = "months"))
nrow(CJ(ids, mons))
# * We have 48 'realised' unit_ids 
# -> These are the individual/months for which we have observations
sort(unique(iteration$unit_id))
length(unique(iteration$unit_id))
# * We have 154 elements in each dataset (one for each _possible_ unit_id)
length(datasets$detections_by_unit)
length(datasets$archival_by_unit)  
length(datasets$behaviour_by_unit)
# * For all AC/ACDC unit_id, we need detections dataset
sapply(unique(iteration$unit_id[iteration$dataset %in% c("ac", "acdc")]), function(id) {
  !is.null(datasets$detections_by_unit[[id]])
}) |> all() |> stopifnot()
# * For all DC/ACDC unit_id, we need archival dataset
sapply(unique(iteration$unit_id[iteration$dataset %in% c("dc", "acdc")]), function(id) {
  !is.null(datasets$archival_by_unit[[id]])
}) |> all() |> stopifnot()
# * For all AC/DC/ACDC unit_id, we need behaviour dataset
sapply(unique(iteration$unit_id), function(id) {
  !is.null(datasets$behaviour_by_unit[[id]])
}) |> all() |> stopifnot()
# Check folder_xinit
all(file.exists(file.path(iteration$folder_xinit, "xinit-fwd.qs"))) |> stopifnot()
all(file.exists(file.path(iteration$folder_xinit, "xinit-bwd.qs"))) |> stopifnot()
# Check number of rows
table(iteration$dataset)
table(iteration$dataset, iteration$mobility)

#### Estimate coordinates
if (FALSE) {
  
  #### Select batch
  # As for the simulations, we batch by mobility
  batch     <- pars$pmovement$mobility[1]
  iteration <- iteration[mobility == batch, ]
  set_vmap(.vmap = here_data("spatial", glue("vmap-{batch}.tif")))
  
  #### Select iterations
  # Select iterations 
  # * Select one AC/DC/ACDC run for testing 
  nrow(iteration)
  # iteration <- 
  #   iteration |> 
  #   filter(sensitivity == "best") |>
  #   group_by(dataset) |> 
  #   slice(1L) |>
  #   as.data.table()
  # * Select by dataset/sensitivity
  iteration <- iteration[dataset == "ac", ]
  # iteration <- iteration[sensitivity == "best", ]
  # * Select rows
  # if (batch == pars$pmovement$mobility[1]) {
  #   rows <- 1:30      # 1
  #   # rows <- 31:60     # 2
  #   # rows <- 61:90     # 3
  #   # rows <- 91:120    # 4
  #   # rows <- 121:150   # 5
  #   # rows <- 151:180   # 6
  #   # rows <- 181:210   # 7
  #   # rows <- 211:240   # 8
  # } else {
  #   rows <- 1:48      # 9
  #   rows <- 1:48      # 10
  # }
  rows <- 1:24
  # rows <- 45:48
  iteration <- iteration[rows, ]
  
  #### Estimate coordinates
  # * Convergence trials complete (see /log/real/trials/log-summary.txt)
  # * AC:
  # - Batch 1 (144), SIA-USER024-P, 8 threads, 3.6 days
  # - Batch 2 (48), siam-linux20
  # - Batch 3 (48), siam-linux20
  # * DC:
  # - Batch 1 (144): MCC02XT0AZJGH5, 12 threads, 41 hr; 10 threads, 42.9 hour (incl. 15 min break every 2 hr)
  # - Batch 2 (48): siam-linux20
  # - Batch 3 (48): siam-linux20
  # * ACDC:
  # - Batch 1 (240): siam-linux20
  # - Batch 2 (48): siam-linux20
  # - Batch 3 (48): siam-linux20
  gc()
  nrow(iteration)
  lapply_estimate_coord_patter(iteration  = iteration,
                               datasets   = datasets, 
                               trial      = FALSE, 
                               log.folder = here_data("output", "log", "real"), 
                               log.txt    = glue("log-{iteration$dataset[1]}-{batch}-{min(rows)}.txt"))
  
}

#### Patter error check 
# (Code copied from simulate-algorithms.R)
# Check for errors on forward filter 
sapply(split(iteration, seq_row(iteration)), function(d) {
  qs::qread(file.path(d$folder_coord, "data-fwd.qs"))$error
}) |> unlist() |> unique()
# Check for errors on backward filter 
sapply(split(iteration, seq_row(iteration)), function(d) {
  file <- file.path(d$folder_coord, "data-bwd.qs")
  if (file.exists(file)) {
    qs::qread(file)$error
  }
}) |> unlist() |> unique()
# Check for errors on smoother
sapply(split(iteration, seq_row(iteration)), function(d) {
  file <- file.path(d$folder_coord, "data-smo.qs")
  if (file.exists(file)) {
    qs::qread(file)$error
  }
}) |> unlist() |> unique()

#### Estimate UDs
if (FALSE && patter:::os_linux()) {
  
  # Examine selected coord
  # lapply_qplot_coord(iteration, 
  #                    "coord-smo.qs",
  #                    extract_coord = function(s) s$states[sample.int(1000, size = .N, replace = TRUE), ])
  
  #### Estimate UDs
  # Time trial 
  lapply_estimate_ud_spatstat(iteration     = iteration[1, ], 
                              extract_coord = function(s) s$states,
                              cl            = NULL, 
                              plot          = FALSE)
  # Implementation (~43 min)
  lapply_estimate_ud_spatstat(iteration     = iteration, 
                              extract_coord = function(s) s$states,
                              cl            = 10L, 
                              plot          = FALSE)
  # (optional) Examine selected UDs
  lapply_qplot_ud(iteration, "spatstat", "h", "ud.tif")
  
}


###########################
###########################
#### Generate outputs


###########################
#### Tidy iteration (copied from simulate-algorithms.R)

# Define grouping variables with tidy labels
iteration[, unit_id := factor(unit_id)]
iteration[, dataset := factor(dataset, levels = c("ac", "dc", "acdc"), labels = c("AC", "DC", "ACDC"))]
# Revise 'sensitivity' coding
iteration[sensitivity == "movement" & mobility == pars$pmovement$mobility[2], sensitivity := "movement-under"]
iteration[sensitivity == "movement" & mobility == pars$pmovement$mobility[3], sensitivity := "movement-over"]
iteration[sensitivity == "ac" & receiver_alpha == pars$pdetection$receiver_alpha[2], sensitivity := "ac-under"]
iteration[sensitivity == "ac" & receiver_alpha == pars$pdetection$receiver_alpha[3], sensitivity := "ac-over"]
iteration[sensitivity == "dc" & depth_sigma == pars$pdepth$depth_sigma[2], sensitivity := "dc-under"]
iteration[sensitivity == "dc" & depth_sigma == pars$pdepth$depth_sigma[3], sensitivity := "dc-over"]
iteration[, sensitivity := factor(sensitivity, 
                                  levels = c("best", 
                                             "movement-under", "movement-over", 
                                             "ac-under", "ac-over",
                                             "dc-under", "dc-over"), 
                                  labels = c("Best", 
                                             "Move (-)", "Move (+)", 
                                             "AC(-)", "AC(+)", 
                                             "DC(-)", "DC(+)"))] 


###########################
#### Convergence 

# Determine convergence
iteration[, convergence := sapply(split(iteration, seq_row(iteration)), function(d) {
  file.exists(file.path(d$folder_coord, "coord-smo.qs"))
})]

# Summarise convergence
iteration |> 
  group_by(dataset, sensitivity) |> 
  summarise(convergence = length(which(convergence)) / n())

# Plot convergence proportions 
iteration |> 
  group_by(dataset, sensitivity) |> 
  summarise(convergence = length(which(convergence)) / n()) |>
  ungroup() |>
  as.data.table() |> 
  ggplot(aes(dataset, convergence, fill = sensitivity)) +
  geom_bar(stat = "identity", 
           position = "dodge",
           colour = "black", linewidth = 0.2, 
           alpha = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
  labs(
    x = "Algorithm",
    y = "Pr(convergence)",
    fill = "Algorithm"
  ) + 
  theme_bw() + 
  theme(axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        panel.grid = element_blank())


###########################
#### Mapping (~1 min)

if (FALSE) {
  
  unique(iteration$dataset)
  unique(iteration$sensitivity)
  length(unique(iteration$individual_id))
  length(unique(iteration$month_id))
  length(unique(paste(iteration$sensitivity, iteration$dataset)))
  cl_lapply(unique(iteration$sensitivity), .cl = 9L, .fun = function(sens){
    lapply(c("AC", "DC", "ACDC"), function(alg) {
      
      # Define mapfiles
      mapfiles <-
        iteration |> 
        filter(sensitivity == sens & dataset == alg) |>
        mutate(mapfile = file.path(folder_ud, "spatstat", "h", "ud.tif"), 
               individual_id = factor(individual_id, levels = sort(unique(individual_id)))) |>
        dplyr::select(row = individual_id, column = month_id, mapfile) |>
        as.data.table() 
      
      # Make map
      if (nrow(mapfiles)) {
        png_args <- 
          list(filename = here_fig("analysis", glue("map-patter-{alg}-{sens}.png")), 
               height = 10, width = 9, units = "in", res = 800)
        ggplot_maps(mapdt = mapfiles, png_args = png_args)
      }
      NULL
    })
  })
  
}


###########################
#### Compute residency (~22 s)

iteration_res <- copy(iteration)
iteration_res[, file := file.path(iteration$folder_ud, "spatstat", "h", "ud.tif")]
iteration_res[, file_exists := file.exists(file)]
residency <- lapply_estimate_residency_ud(files = iteration_res$file[iteration_res$file_exists])
# Write output
residency <- 
  left_join(iteration_res, residency, by = "file") |> 
  mutate(algorithm = dataset) |>
  select(individual_id, month_id, unit_id, algorithm, sensitivity, zone, time) |> 
  arrange(individual_id, month_id, unit_id, algorithm, sensitivity, zone) |>
  as.data.table()
qs::qsave(residency, here_data("output", "analysis-summary", "residency-patter.qs"))


#### End of code.
###########################
###########################