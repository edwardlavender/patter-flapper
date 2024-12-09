if (!patter:::os_linux() | Sys.getenv("JULIA_SESSION") == "FALSE") {

#' Set par parameters
set_par <- function(...) {
  par(mgp = c(3, 0.4, 0), tcl = -0.25, ...)
}

#' Add distribution polygons
add_dbn <- function(x, y,
                    col = scales::alpha("lightgrey", 0.25), ...) {
  polygon(c(x, rev(x)), c(y, rep(0, length(y))), col = col, ...)
  invisible(NULL)
}

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

}
