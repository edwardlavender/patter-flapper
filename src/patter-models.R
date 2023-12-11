#' @title Depth error function
# https://github.com/edwardlavender/flapper_appl/blob/master/R/define_global_param.R

calc_depth_error <- readRDS(dv::here_data("input", "pars.rds"))$patter$calc_depth_error

#' @title Origin spatRaster for *DC models (DCPF, ACDCPF)

dc_origin <- function(.bathy, .depth, .calc_depth_error) {
  error   <- .calc_depth_error(.depth)
  shallow <- .depth - (error + 12.5) 
  deep    <- .depth + error
  terra::clamp(.bathy, shallow, deep, values = FALSE)
}

#' @title Depth weighting function

update_ac <- function(.particles, .bathy, .obs, .t, ...) {
  # Extract depth of seabed at particle positions
  if (!rlang::has_name(.particles, "bathy")) {
    .particles[, bathy := terra::extract(.bathy, .particles$cell_now)]
  }
  # Define bathymetry data source (digi versus howe) & define adjustment for depth-error model
  # > For some receivers in the Howe data, the individual is ~10-20 m shallower than previously assumed.
  # > For the Digimap data, the individual can be up to 40 m shallower. 
  .particles[, digi := terra::extract(.bset, .particles$cell_now)]
  .particles[, boost := 20L]
  bool_digi <- .particles$digi == 1L
  if (any(.particles$digi == 1L)) {
    .particles[which(bool_digi), boost := 40L]
  }
  # Define shallow and deep depth limits for each particle based on location
  .particles[, error := calc_depth_error(depth)]
  .particles[, deep := depth + error]
  .particles[, shallow_1 := depth - (error + 5)]
  .particles[, shallow_2 := depth - (error + boost)]
  # Define particle weights
  weight <- rep(0L, nrow(.particles))
  pos <- which(.particles$bathy <= .particles$deep & 
                 .particles$bathy  >= .particles$shallow_1)
  weight[pos] <- 1L
  pos <- which(.particles$bathy < .particles$shallow_1 & 
                 .particles$bathy >= .particles$shallow_2)
  weight[pos] <- 0.01
  weight / sum(weight)
  
}

#' @title Plot the depth-error model

add_depth_error_model <- function(depth) {
  # Define depth limits
  error <- calc_depth_error(depth)
  deep <- depth + error
  boost_howe <- 20
  boost_digi <- 40
  shallow_1 <- depth - (error + 5)
  shallow_2_howe <- depth - (error + boost_howe)
  shallow_2_digi <- depth - (error + boost_digi)
  # Scale depth limits
  if (TRUE) {
    deep <- deep - depth
    shallow_1 <- shallow_1 - depth 
    shallow_2_howe <- shallow_2_howe - depth
    shallow_2_digi <- shallow_2_digi - depth
  }
  # Add lines to a plot
  lines(c(shallow_2_digi, shallow_2_digi), c(0, 0.01), col = "red")
  lines(c(shallow_2_howe, shallow_2_howe), c(0, 0.01))
  lines(c(shallow_2_digi, shallow_2_howe), c(0.01, 0.01), col = "red")
  lines(c(shallow_2_howe, shallow_1), c(0.01, 0.01))
  lines(c(shallow_1, shallow_1), c(0.01, 1))
  lines(c(shallow_1, deep), c(1, 1))
  lines(c(deep, deep), c(1, 0))
}
