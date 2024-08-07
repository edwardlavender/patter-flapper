#' Plot selected coords

lapply_qplot_coord <- function(iteration, ..., extract_coord = NULL) {
  bathy <- terra::rast(here_data("spatial", "bathy.tif"))
  pp <- par(mfrow = c(2, 2))
  sapply(1:4, function(i) {
    try({
      print(i)
      coord.qs <- file.path(iteration$folder_coord[i], ...)
      coord    <- qs::qread(coord.qs)
      if (!is.null(extract_coord)) {
        coord <- extract_coord(coord)
      }
      terra::plot(bathy, maxcell = 100000, main = i)
      points(coord$x, coord$y, pch = ".", col = "red")
    })
    invisible(NULL)
  })
  par(pp)
  invisible(NULL)
}

#' Plot selected UDs
lapply_qplot_ud <- function(iteration, ...) {
  pp <- par(mfrow = c(2, 2))
  sapply(1:4, function(i) {
    try({
      print(i)
      ud.tif <- file.path(iteration$folder_ud[i], ...)
      ud     <- terra::rast(ud.tif)
      terra::plot(ud, main = i)
    })
    invisible(NULL)
  })
  par(pp)
  invisible(NULL)
}