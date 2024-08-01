#' @title qread/save helpers for terra

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

#' @title Assign SpatRasters in Julia

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
  JuliaCall::julia_command(glue::glue('{x} = GeoArrays.read("{file}");'))
  invisible(NULL)
}

