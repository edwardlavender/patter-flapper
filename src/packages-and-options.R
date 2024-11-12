# Packages
library(collapse)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(dv)
library(ggplot2)
library(JuliaCall)
library(lubridate)
library(patter)
library(spatstat)
library(tictoc)

# terra options
if (!patter:::os_linux()) {
  op <- options(terra.pal = rev(terrain.colors(256)))
  if (Sys.info()["nodename"] == "MCC02XT0AZJGH5") {
    tmpdir <- "/Volumes/My Book/projects/eawag/patter-flapper/data-raw/temporary"
    stopifnot(dir.exists(tmpdir))
    top <- terra::terraOptions(tempdir = tmpdir)
  }
}

# spatstat options
# * Use npixel to set the pixel resolution 
# * It is essential library(spatstat) is called above for this to work properly
# * For testing, set npixel to small (e.g., 50 pixel)
# * Otherwise, finer resolution is better
# * The pixel resolution should just be high enough for the plot to look nice &
# ... low enough for the calculation not to take too long.
# See https://stackoverflow.com/questions/67687984/choosing-pixel-size-for-kernel-density-estimate-map-density-ppp-in-r 
sop <- spatstat.geom::spatstat.options(npixel = 500)
sop$npixel

# Multi-threading
data.table::setDTthreads(threads = Sys.getenv("OMP_NUM_THREADS"))
data.table::getDTthreads()

# Progress bar (chunk) options for cl_lapply()
# * This is the number of chunks per core
# * We read rasters into memory on each chunk (slow), so this number is relatively low
pbo <- pbapply::pboptions("nout" = 10)
