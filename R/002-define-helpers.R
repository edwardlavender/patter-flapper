###########################
###########################
#### define_helpers.R

#### Aims
# 1) Define helper functions

#### Prerequisites
# 1) NA


###########################
###########################
#### Define functions

#' @title Depth error function
# https://github.com/edwardlavender/flapper_appl/blob/master/R/define_global_param.R

calc_depth_error <- function(depth) {
  e <- 4.77 + 2.5 + sqrt(0.5^2 + (0.013 * depth)^2)
  matrix(c(-(e + 5), e), nrow = 2)
}
calc_depth_error <- Vectorize(calc_depth_error)

#' @title Depth weighting function

update_ac <- function(.particles, .bathy, .obs, .t, ...) {
  # Extract depth of seabed at particle positions
  .particles$bathy <- terra::extract(.bathy, as.matrix(.particles[, c("x_now", "y_now")]))
  # Weight = 1 in locations where bathymetric depth is within possible limits, otherwise 0
  (.particles$bathy  >= .obs$depth_shallow[.t] & .particles$bathy <= .obs$depth_deep[.t]) + 0
}


#### End of code.
###########################
###########################