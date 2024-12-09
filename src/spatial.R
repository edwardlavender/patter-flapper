if (!patter:::os_linux() | Sys.getenv("JULIA_SESSION") == "FALSE") {

#' qread/qsave helpers for terra
qreadvect <- function(...) {
  qs::qread(...) |> terra::vect()
}
qsavevect <- function(x, file, ...) {
  x |> 
    sf::st_as_sf() |> 
    qs::qsave(file = file, ...)
}
qreadext <- function(...) {
  qs::qread(...) |> terra::ext()
}
qsaveext <- function(x, file, ...) {
  c(x[1], x[2], x[3], x[4]) |> 
    qs::qsave(file = file, ...)
}

#' Assign SpatRasters in Julia
julia_assign_SpatRaster <- function(x, value) {
  # Check SpatRaster
  stopifnot(inherits(value, "SpatRaster"))
  # Define file
  file <- terra::sources(value)
  if (file == "") {
    file <- tempfile(fileext = ".tif")
    terra::writeRaster(value, file)
  }
  # Set env
  file <- normalizePath(file, winslash = "/", mustWork = TRUE)
  JuliaCall::julia_command(glue::glue('{x} = GeoArrays.read("{file}");'))
  invisible(NULL)
}

#' Normalise SpatRasters
spatNormalise <- function(x) {
  x / as.numeric(terra::global(x, "sum", na.rm = TRUE)[, 1])
}

#' Distance functions
#' https://github.com/edwardlavender/patter/blob/70e69537ea71e5a77b18e5457d6c1a0beb68190a/R/dists.R 
dist_along_path <- function(.xy, .lonlat = FALSE) {
  if (!inherits(.xy, "matrix")) {
    if (ncol(.xy) != 2L) {
      abort("`.xy` should be a two-column matrix of coordinates.")
    }
    .xy <- as.matrix(.xy, ncol = 2)
  }
  dist <- terra::distance(.xy,
                          lonlat = .lonlat,
                          sequential = TRUE)
  c(dist[2:length(dist)], NA_real_)
}

}
