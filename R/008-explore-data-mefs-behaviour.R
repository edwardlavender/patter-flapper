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
#### Analyse resting behaviour

#### Define 'resting' versus 'active' states, as in process-data-mefs
archival[, fct := paste(individual_id, mmyy)]
archival[, state := get_mvt_resting(copy(archival), fct = "fct")]

#### Calculate the duration of each state 
# Duration of resting behaviour 
resting <- duration_state(archival, fct = "fct", state = 0L)
# Duration of active behaviour 
active <- duration_state(archival, fct = "fct", state = 1L)
# > NB: Duration here is the number of two-min time steps for which the behaviour was maintained

#### Summarises
# Basic stats on the proportions of time spent resting versus active by individual
archival |> 
  group_by(individual_id) |> 
  summarise(pr = length(which(state == 0)) / n()) |> 
  pull(pr) |> 
  utils.add::basic_stats()
# Time spent resting 
hist(resting$duration)
hist(active$duration)

#### Model the duration of resting or active behaviour

# Define data
data <- resting$duration
# data <- active$duration

# Explore suitable distributions
fitdistrplus::descdist(data)

# Fit a distribution to the data
# > None of these distributions are a perfect fit to the data
# > We will use the gamma distribution
dbn <- "gamma"; ddbn <- "gamma"
# dbn <- "log-normal"; ddbn <- "lnorm" 
# dbn <- "weibull"; ddbn <- "weibull"
# dbn <- "poisson"; ddbn <- "pois"
# dbn <- "negative binomial"; ddbn <- "nbinom"
(pars <- MASS::fitdistr(data, dbn))
fitdistrplus::plotdist(data, ddbn, as.list(pars$estimate))

# Resting parameters
# shape           rate    
# 0.7297834831   0.1626763493 
# (0.0026985241) (0.0008370212)

# Active parameters
# shape           rate    
# 0.9603830324   0.1679107113 
# (0.0036422085) (0.0008244781)

# hist(data, breaks = 5, probability = TRUE)
# curve(dgamma(x, pars$estimate[1], pars$estimate[2]), col = "red", lwd = 2, add = TRUE)

# Simulate states
timeline <- seq(as.POSIXct("2024-04-01 00:00:00", tz = "UTC"), 
                as.POSIXct("2024-04-30 23:58:00", tz = "UTC"), 
                by = "2 mins")
simulate_behaviour(timeline)


#### End of code. 
###########################
###########################