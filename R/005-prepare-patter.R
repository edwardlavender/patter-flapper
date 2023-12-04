###########################
###########################
#### prepare-patter.R

#### Aims
# 1) Prepare flapper algorithm inputs

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

#### Detection overlaps (~0 s)
tic()
overlaps <- acs_setup_detection_overlaps(moorings)
toc()

#### Detection kernels
# * We can use the default settings here. 
# * Note that MEFS::moorings contains 48 receivers
# * Whereas patter::dat_moorings only contains a subset of those receivers
# Timings:
# * receiver-specific inverse kernels: ~90 mins 
# * area-wide kernels: ~24 min
# * total: 115 mins
tic()
kernels <- acs_setup_detection_kernels(moorings, 
                                       .calc_detection_pr = acs_setup_detection_pr, 
                                       .bathy = bathy)
toc()

### Save outputs
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


###########################
###########################
#### Prepare {spatstat} inputs

#### Define image (~8 s)
tic()
im <- as.im.SpatRaster(bathy)
toc()

#### Define observation window (~2 s + 164 s)
tic()
win <- as.owin.sf(coast, .invert = TRUE)
toc()
tic()
plot(win, col = "blue")
toc()

#### Save outputs (~6 s + 115 s)
tic()
qs::qsave(im, here_data("spatial", "im.qs"))
toc()
tic()
qs::qsave(win, here_data("spatial", "win.qs"))
toc()


#### End of code. 
###########################
###########################