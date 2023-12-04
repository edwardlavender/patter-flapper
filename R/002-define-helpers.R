###########################
###########################
#### define_helpers.R

#### Aims
# 1) Define helper functions

#### Prerequisites
# 1) NA


###########################
###########################
#### Utilities

#' @title {Tools4ETS} imports

difference <- function(x2, x1, f = NULL, ...) {
  if (class(x2)[1] %in% c("numeric", "integer")) {
    d <- x2 - x1
  }
  else if (class(x2)[1] %in% c("POSIXct", "POSIXlt", "Date")) {
    d <- difftime(x2, x1, ...)
  }
  if (!is.null(f)) {
    d <- as.numeric(d)
  }
  d
}

serial_difference <- function(x, na.rm = FALSE, ...) {
  dur <- difference(dplyr::lead(x), x, ...)
  if (na.rm) {
    posNA <- which(is.na(dur))
    dur <- dur[-c(posNA)]
  }
  dur
}


###########################
###########################
#### Spatial helpers

#' @title Spatial helpers

readRasterLs <- function(folder, ...) {
  files <- list.files(folder, full.names = TRUE)
  ids   <- tools::file_path_sans_ext(basename(files))
  rasters <- pbapply::pblapply(files, function(file) {
    terra::rast(file, ...)
  })
  names(rasters) <- ids
  index <- as.character(seq_len(max(as.integer(ids))))
  lapply(index, function(i) {
    rasters[[i]]
  })
}

writeRasterLs <- function(x, folder, ...) {
  index <- seq_len(length(x))
  outfiles <- file.path(folder, paste0(index, ".tif"))
  pbapply::pblapply(index, function(i) {
    if (!is.null(x[[i]])) {
      terra::writeRaster(x[[i]], outfiles[i], ...)
    }
  })
  invisible(NULL)
}

###########################
###########################
#### {patter} helpers

#' @title {patter} helpers

as_glatos <- function(.acoustics) {
  data.frame(detection_timestamp_utc = .acoustics$timestamp, 
             transmitter_codespace = "000",
             transmitter_id = as.character(.acoustics$individual_id), 
             receiver_sn = as.character(.acoustics$receiver_id)
  )
}

acs_setup_detection_kernels_read <- function() {
  kernels <- list()
  kernels$array_design <- 
    readRDS(here_data("input", "kernels", "array_design.rds"))
  kernels$array_design_by_date <- 
    readRDS(here_data("input", "kernels", "array_design_by_date.rds"))
  kernels$receiver_specific_kernels <- 
    readRasterLs(here_data("input", "kernels", "receiver-specific-kernels"))
  kernels$receiver_specific_inv_kernels <- 
    readRasterLs(here_data("input", "kernels", "receiver-specific-inv-kernels"))
  kernels$bkg_surface_by_design <- 
    readRasterLs(here_data("input", "kernels", "bkg-surface-by-design"))
  kernels$bkg_inv_surface_by_design <- 
    readRasterLs(here_data("input", "kernels", "bkg-inv-surface-by-design"))
  kernels
}

# (optional) TO DO
# * Add acs_setup_detection_kernels_write()

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
  if (!rlang::has_name(.particles, "bathy")) {
    .particles[, bathy := terra::extract(.bathy, .particles$cell_now)]
  }
  # Weight = 1 in locations where bathymetric depth is within possible limits, otherwise 0
  (.particles$bathy  >= .obs$depth_shallow[.t] & .particles$bathy <= .obs$depth_deep[.t]) + 0
}


#### End of code.
###########################
###########################