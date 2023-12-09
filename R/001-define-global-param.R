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


###########################
###########################
#### Define parameters

# Detection range
detection_range <- 750

# Mobility 
mobility <- 500

# (original) Depth-error function 
calc_depth_error_1 <- function(depth) {
  e <- 4.77 + 2.5 + sqrt(0.5^2 + (0.013 * depth)^2)
  matrix(c(-(e + 5), e), nrow = 2)
}
calc_depth_error_1 <- Vectorize(calc_depth_error_1)

# (revised) Depth-error function 
calc_depth_error_2 <- function(depth) {
  4.77 + 2.5 + sqrt(0.5^2 + (0.013 * depth)^2)
}

# Collate pars
pars <- 
  list(
    flapper = list(
      detection_range = detection_range, 
      mobility = mobility, 
      calc_depth_error = calc_depth_error_1), 
    patter = list(
      detection_range = detection_range, 
      mobility = mobility, 
      calc_depth_error = calc_depth_error_2)
  )

# Save pars
saveRDS(pars, here_data("input", "pars.rds"))


#### End of code. 
###########################
###########################