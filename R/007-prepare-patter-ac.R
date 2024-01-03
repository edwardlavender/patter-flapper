###########################
###########################
#### prepare-patter-ac.R

#### Aims
# 1) Prepare *ACPF algorithm components
#    * Detection overlaps
#    * Detection kernels

#### Prerequisites
# 1) Obtain raw data
# 2) https://github.com/edwardlavender/flapper_appl/blob/master/R/define_global_param.R


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
dv::src()
library(dv)
library(patter)
library(tictoc)

#### Load data
coast    <- readRDS(here_data("spatial", "coast.rds"))
bathy    <- terra::rast(here_data("spatial", "bathy.tif"))
moorings <- readRDS(here_data("mefs", "moorings.rds"))

#### Local variables
overwrite <- FALSE


###########################
###########################
#### Prepare AC* inputs

#### Data list
dlist <- pat_setup_data(.moorings = moorings, .bathy = bathy, .lonlat = FALSE)

#### Detection overlaps (~0 s)
tic()
overlaps <- acs_setup_detection_overlaps(dlist)
toc()

#### Detection kernels
# * We can use the default settings here. 
# * Note that MEFS::moorings contains 48 receivers
# * Whereas patter::dat_moorings only contains a subset of those receivers
# Timings:
# * receiver-specific kernels: ~3 s
# * area-wide kernels:  ~51 s
# * processing (extension & masking): 
# * total: ~42 mins
tic()
kernels <- acs_setup_detection_kernels(dlist)
toc()

#### Save outputs
# ~20 mins (~7 mins per raster list)
tic()
saveRDS(overlaps, here_data("input", "overlaps.rds"))
saveRDS(kernels$array_design, 
        here_data("input", "kernels", "array_design.rds"))
saveRDS(kernels$array_design_by_date, 
        here_data("input", "kernels", "array_design_by_date.rds"))
writeRasterLs(kernels$receiver_specific_kernels,
              here_data("input", "kernels", "receiver-specific-kernels"))
writeRasterLs(kernels$receiver_specific_inv_kernels,
              here_data("input", "kernels", "receiver-specific-inv-kernels"))
writeRasterLs(kernels$bkg_surface_by_design,
              here_data("input", "kernels", "bkg-surface-by-design"))
writeRasterLs(kernels$bkg_inv_surface_by_design,
              here_data("input", "kernels", "bkg-inv-surface-by-design"))
toc()

#### Validate revised acs_setup_detection_kernels() routine (revised on 2024-01-03) (~0 s)
# kernels <- acs_setup_detection_kernels_read()
old         <- here_data("input", "kernels", "2024-01-03")
old_kernels <- acs_setup_detection_kernels_read(path = old)
all.equal(kernels$array_design, old_kernels$array_design)
all.equal(kernels$array_design_by_date, old_kernels$array_design_by_date)
spatEqual <- function(x, y) {
  # print(x)
  if (!is.null(x)) {
    terra::compareGeom(x, y, stopOnError = TRUE)
    terra::ext(x) <- terra::ext(bathy)
    terra::ext(y) <- terra::ext(bathy)
    stopifnot(terra::all.equal(x, y, maxcell = 1e5L, tolerance = 1e-6))
    # Too slow: 
    # stopifnot(all(abs(x - y) < 1e-6))
  }
  invisible()
}
# The revised approach returns identical outputs, up to a tolerance of 1e-6
# * The origin and spatExtent are not (quite) identical
# * But very close
pbapply::pbmapply(spatEqual,
                  kernels$receiver_specific_kernels, 
                  old_kernels$receiver_specific_kernels) |> invisible()
pbapply::pbmapply(spatEqual, 
                  kernels$receiver_specific_inv_kernels, 
                  old_kernels$receiver_specific_inv_kernels) |> invisible()
pbapply::pbmapply(spatEqual, 
                  kernels$bkg_surface_by_design, 
                  old_kernels$bkg_surface_by_design) |> invisible()
pbapply::pbmapply(spatEqual, 
                  kernels$bkg_inv_surface_by_design, 
                  old_kernels$bkg_inv_surface_by_design) |> invisible()


#### End of code. 
###########################
###########################