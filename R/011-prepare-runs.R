###########################
###########################
#### prepare-runs.R

#### Aims
# 1) Prepare folders & iteration datasets for algorithm runs

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
unitsets          <- qs::qread(here_data("input", "unitsets.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))
archival_by_unit  <- qs::qread(here_data("input", "archival_by_unit.qs"))
pars              <- qs::qread(here_data("input", "pars.qs"))


###########################
###########################
#### Define base folders

tic()

# (optional) Clean data/output/analysis
if (FALSE) {
  unlink(file.path("data", "output", "analysis"), recursive = TRUE)
}

# {individual}/{month}/{package}/
unitsets[, folder_unit := file.path("data", "output", "analysis", 
                                     individual_id, month_id)]
dirs.create(unitsets$folder_unit)


###########################
###########################
#### Define COA iteration dataset

#### Identify individuals/months with detections
iteration_ac <- 
  unitsets |> 
  filter(dataset == "ac") |> 
  as.data.table()

#### Define COA parameters (based on simulations)
# > We consider a 'best', 'restrictive' and 'flexible' parameter value
parameters <- data.table(parameter_id = 1:3L, 
                         delta_t = c("2 days", "1 day", "3 days"))

#### Define iteration dataset
iteration_coa <- 
  cross_join(iteration_ac, parameters) |> 
  arrange(unit_id, parameter_id) |> 
  mutate(index = row_number(), 
         folder_coord = file.path(folder_unit, "coa", dataset, parameter_id, "coord"),
         folder_ud = file.path(folder_unit, "coa", dataset, parameter_id, "ud")) |> 
  select("index", 
         "unit_id", "individual_id",  "month_id", 
         "dataset", 
         "parameter_id", "delta_t", "folder_coord", "folder_ud") |> 
  as.data.table()

#### Build COA folders
dirs.create(iteration_coa$folder_coord)
dirs.create(iteration_coa$folder_ud)
dirs.create(file.path(iteration_coa$folder_ud, "spatstat", "h"))


###########################
###########################
#### Define RSP iteration dataset

#### Identify individuals/months with detections
# > Implemented above. 

#### Define parameters (based on simulations)
parameters <- data.table(parameter_id = 1:3L, 
                         er.ad = c(500, 250, 750))

#### Define iteration dataset
iteration_rsp <- 
  cross_join(iteration_ac, parameters) |> 
  arrange(unit_id, parameter_id) |> 
  mutate(index = row_number(), 
         folder_coord = file.path(folder_unit, "rsp", dataset, parameter_id, "coord"), 
         folder_ud = file.path(folder_unit, "rsp", dataset, parameter_id, "ud")) |> 
  select("index", 
         "unit_id", "individual_id",  "month_id", 
         "dataset", 
         "parameter_id", "er.ad", 
         "folder_coord", "folder_ud") |> 
  as.data.table()

#### Build RSP folders
dirs.create(iteration_rsp$folder_coord)
dirs.create(iteration_rsp$folder_ud)
dirs.create(file.path(iteration_rsp$folder_ud, "spatstat", "h"))
dirs.create(file.path(iteration_rsp$folder_ud, "dbbmm"))


###########################
###########################
#### Define patter iteration dataset

#### Define parameters
# Define parameters for each process, including:
# * 'best' guess parameters (first value)
# * 'restrictive' parameters (second value)
# * 'flexible' parameters (third value)
pmovement  <- pars$pmovement
pdetection <- pars$pdetection
pdepth     <- pars$pdepth
# Collate parameter combinations
# * Best guess parameters
# * Variable movement parameters (with others held at best-guess values)
# * Variable detection parameters (with others held at best-guess values)
# * Variable depth parameters (with others held at best-guess values)
parameters <- 
  rbind(
    cbind(sensitivity = "best", pmovement[1, ], pdetection[1, ], pdepth[1, ]),
    cbind(sensitivity = "movement", pmovement[2:3, ], pdetection[1, ], pdepth[1, ]),
    cbind(sensitivity = "ac", pmovement[1, ], pdetection[2:3, ], pdepth[1, ]),
    cbind(sensitivity = "dc", pmovement[1, ], pdetection[1, ], pdepth[2:3, ]) 
  )
parameters[, parameter_id := seq_len(.N)]

#### Define iteration dataset
iteration_patter <- lapply(split(unitsets, seq_len(nrow(unitsets))), function(d) {
  
  # Keep the relevant parameters, dependent upon the algorithm
  if (d$dataset == "ac") {
    p <- parameters[sensitivity %in% c("best", "movement", "ac"), ]
    p[, c("depth_sigma", "depth_deep_eps") := NA_real_]
  }
  if (d$dataset == "dc") {
    p <- parameters[sensitivity %in% c("best", "movement", "dc"), ]
    p[, c("receiver_alpha", "receiver_beta", "receiver_gamma") := NA_real_]
  }
  if (d$dataset == "acdc") {
    p <- copy(parameters)
  }
  cbind(d, p)
  
}) |> 
  rbindlist() |> 
  mutate(index = row_number(), 
         folder_coord = file.path(folder_unit, "patter", dataset, parameter_id, "coord"), 
         folder_ud = file.path(folder_unit, "patter", dataset, parameter_id, "ud")) |> 
  select("index", 
         "unit_id", "individual_id",  "month_id", 
         "dataset", 
         "parameter_id", "sensitivity", 
         "k1", "k2", "theta1", "theta2", "mobility", 
         "receiver_alpha", "receiver_beta", "receiver_gamma", 
         "depth_sigma", "depth_deep_eps",
         "folder_coord", "folder_ud") |>
  as.data.table()
# Add additional controls
np <- data.table(dataset = c("ac", "dc", "acdc"), 
                 np = c(1.5e5L, 1.5e5L, 1.5e5L))
iteration_patter[, np := np$np[match(dataset, np$dataset)]]
iteration_patter[, smooth := TRUE]

#### Build patter folders
dirs.create(iteration_patter$folder_coord)
dirs.create(iteration_patter$folder_ud)
dirs.create(file.path(iteration_patter$folder_ud, "spatstat", "h"))
dirs.create(file.path(iteration_patter$folder_ud, "pou"))


###########################
###########################
#### Write outputs to file

qs::qsave(iteration_coa, here_data("input", "iteration", "coa.qs"))
qs::qsave(iteration_rsp, here_data("input", "iteration", "rsp.qs"))
qs::qsave(iteration_patter, here_data("input", "iteration", "patter.qs"))


toc()


#### End of code.
###########################
###########################