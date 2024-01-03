#' @title {patter} helpers

# Convert .acoustics to glatos format
as_glatos <- function(.acoustics) {
  data.frame(detection_timestamp_utc = .acoustics$timestamp, 
             transmitter_codespace = "000",
             transmitter_id = as.character(.acoustics$individual_id), 
             receiver_sn = as.character(.acoustics$receiver_id)
  )
}

# Read detection kernels
acs_setup_detection_kernels_read <- function(path = here_data("input", "kernels")) {
  kernels <- list()
  kernels$array_design <- 
    readRDS(file.path(path, "array_design.rds"))
  kernels$array_design_by_date <- 
    readRDS(file.path(path, "array_design_by_date.rds"))
  kernels$receiver_specific_kernels <- 
    readRasterLs(file.path(path, "receiver-specific-kernels"))
  kernels$receiver_specific_inv_kernels <- 
    readRasterLs(file.path(path, "receiver-specific-inv-kernels"))
  kernels$bkg_surface_by_design <- 
    readRasterLs(file.path(path, "bkg-surface-by-design"))
  kernels$bkg_inv_surface_by_design <- 
    readRasterLs(file.path(path, "bkg-inv-surface-by-design"))
  kernels
}

# (optional) TO DO
# * Add acs_setup_detection_kernels_write()