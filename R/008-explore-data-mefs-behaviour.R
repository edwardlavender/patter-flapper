###########################
###########################
#### explore-data-mefs-resting.R

#### Aims
# 1) Examine MEFS datasets (resting behaviour)

#### Prerequisites
# 1) Process MEFS data


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(flapper)
dv::src()

#### Load data
archival  <- qs::qread(here_data("input", "archival_by_unit.qs")) |> plyr::compact() |> rbindlist()


###########################
###########################
#### 

#### Define 'resting' versus 'active' states
archival <- 
  archival <- rbindlist(archival)
archival[, fct := paste(individual_id, mmyy)]

archival[, state := get_mvt_resting(copy(archival), fct = "fct")]

#### Calculate the duration of each state 
archival |> 
  group

duration <- rle(arc$state)
duration <- duration$lengths[duration$values == 0L]

fitdistrplus::descdist(duration)

dbn <- "gamma"
pars <- MASS::fitdistr(duration, dbn)
plotdist(duration, dbn, as.list(pars$estimate))

# hist(duration, breaks = 5, probability = TRUE)
# curve(dcauchy(x, pars$estimate[1], pars$estimate[2]), col = "red", lwd = 2, add = TRUE)


#### End of code. 
###########################
###########################