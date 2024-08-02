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
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(dv)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
dv::src()

#### Load data 
unitsets          <- qs::qread(here_data("input", "unitsets.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))
archival_by_unit  <- qs::qread(here_data("input", "archival_by_unit.qs"))


###########################
###########################
#### Define base folders

# {individual}/{month}/{package}/
unitsets[, folder_base := file.path("data", "output", "analysis", 
                                     individual_id, month_id)]
dirs.create(unitsets$folder_base)


###########################
###########################
#### Define COA iteration dataset

#### Identify individuals/months with detections
iteration_ac <- 
  unitsets |> 
  filter(dataset == "ac") |> 
  as.data.table()

#### Define COA parameters
# > We consider a 'best', 'restrictive' and 'flexible' parameter value
parameters <- data.table(parameter_id = 1:3L, 
                         delta_t = c("6 hours", "12 hours", "24 hours"))

#### Define iteration dataset
iteration_coa <- 
  cross_join(iteration_ac, parameters) |> 
  arrange(unit_id, parameter_id) |> 
  as.data.table()

#### Build COA folders
iteration_coa[, folder := file.path(folder_base, "coa", parameter_id)]
dirs.create(file.path(iteration_coa$folder, "coord"))
dirs.create(file.path(iteration_coa$folder, "ud"))


###########################
###########################
#### Define RSP iteration dataset

#### Identify individuals/months with detections
# > Implemented above. 

#### Define parameters
parameters <- data.table(parameter_id = 1:3L, 
                         er.ad = c(250 * 0.05, 250 * 0.10, 250 * 0.15))

#### Define iteration dataset
iteration_rsp <- 
  cross_join(iteration_ac, parameters) |> 
  arrange(unit_id, parameter_id) |> 
  as.data.table()

#### Build RSP folders
iteration_rsp[, folder := file.path(folder_base, "rsp", parameter_id)]
dirs.create(file.path(iteration_rsp$folder, "coord"))
dirs.create(file.path(iteration_rsp$folder, "ud"))


###########################
###########################
#### Define patter iteration dataset

#### Define parameters
# Define parameters for each process, including:
# * 'best' guess parameters (first value)
# * 'restrictive' parameters (second value)
# * 'flexible' parameters (third value)
pmovement <- data.table(shape = c(15, 15, 15), 
                        scale = c(250, 250, 250), 
                        mobility = c(1000, 500, 1500))
pdetection <- data.table(receiver_alpha = c(4, 5, 2), 
                         receiver_beta = c(-0.01, -0.02, -0.001), 
                         receiver_gamma = c(1500, 1000, 2000))
pdepth <- data.table(depth_sigma = c(50, 10, 100), 
                     depth_deep_eps = c(20, 12, 30))
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
  
}) |> rbindlist()

#### Build patter folders
# {individual}/{mmyy}/{patter}/{algorithm}/{parameter combination}
iteration_patter[, folder := file.path(folder_base, "patter", dataset, parameter_id)]
dirs.create(file.path(iteration_patter$folder, "coord"))
dirs.create(file.path(iteration_patter$folder, "ud"))


###########################
###########################
#### Write outputs to file

qs::qsave(iteration_coa, here_data("input", "iteration", "coa.qs"))
qs::qsave(iteration_rsp, here_data("input", "iteration", "rsp.qs"))
qs::qsave(iteration_patter, here_data("input", "iteration", "patter.qs"))


#### End of code.
###########################
###########################