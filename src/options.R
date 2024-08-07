#' terra colour palette
op <- options(terra.pal = rev(terrain.colors(256)))

#' spatstat options
# For testing, set npixel to small (e.g., 50 pixel)
# Otherwise, finer resolution is better
# The pixel resolution should just be high enough for the plot to look nice &
# ... low enough for the calculation not to take too long.
# See https://stackoverflow.com/questions/67687984/choosing-pixel-size-for-kernel-density-estimate-map-density-ppp-in-r 
sop <- spatstat.geom::spatstat.options(npixel = 50)
sop$npixel

#' Multi-threading
data.table::setDTthreads(threads = Sys.getenv("OMP_NUM_THREADS"))
data.table::getDTthreads()
