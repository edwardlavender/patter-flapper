#' @title Spatial helpers

readRasterLs <- function(folder, ..., index = TRUE) {
  files <- list.files(folder, full.names = TRUE)
  ids   <- tools::file_path_sans_ext(basename(files))
  rasters <- lapply(files, function(file) {
    terra::rast(file, ...)
  })
  names(rasters) <- ids
  if (index) {
    ind <- as.character(seq_len(max(as.integer(ids))))
    rasters <- lapply(ind, function(i) {
      rasters[[i]]
    })
  }
  rasters
}

writeRasterLs <- function(x, folder, ..., index = TRUE) {
  ind <- seq_len(length(x))
  if (index) {
    outfiles <- file.path(folder, paste0(ind, ".tif"))
  } else {
    stopifnot(!is.null(names(x)))
    outfiles <- file.path(folder, paste0(names(x), ".tif"))
  }
  lapply(ind, function(i) {
    if (!is.null(x[[i]])) {
      terra::writeRaster(x[[i]], outfiles[i], ...)
    }
  })
  invisible(NULL)
}
