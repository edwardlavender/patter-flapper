###########################
###########################
#### define-global-param.R

#### Aims
# 1) Define parameters

#### Prerequisites
# 1) flapper_appl project
# 2) These parameters are informed by prepare-patter-param.R


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(dv)


###########################
###########################
#### Define parameters

# step duration
step <- "2 mins"

# Detection range
detection_range <- 750

# Mobility 
mobility <- 500

# (original) Depth-error function 
calc_depth_error <- function(depth) {
  e <- 4.77 + 2.5 + sqrt(0.5^2 + (0.013 * depth)^2)
  matrix(c(-(e + 5), e), nrow = 2)
}
calc_depth_error <- Vectorize(calc_depth_error)

# Collate pars
pars <- 
  list(
    flapper = list(
      step = step,
      detection_range = detection_range, 
      mobility = mobility, 
      calc_depth_error = calc_depth_error), 
    patter = list(
      step = step,
      detection_range = detection_range, 
      mobility = mobility)
  )

# Save pars
saveRDS(pars, here_data("input", "pars.rds"))


#### End of code. 
###########################
###########################