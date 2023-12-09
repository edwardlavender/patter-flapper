#' @title Spatial helpers

readRasterLs <- function(folder, ...) {
  files <- list.files(folder, full.names = TRUE)
  ids   <- tools::file_path_sans_ext(basename(files))
  rasters <- pbapply::pblapply(files, function(file) {
    terra::rast(file, ...)
  })
  names(rasters) <- ids
  index <- as.character(seq_len(max(as.integer(ids))))
  lapply(index, function(i) {
    rasters[[i]]
  })
}

writeRasterLs <- function(x, folder, ...) {
  index <- seq_len(length(x))
  outfiles <- file.path(folder, paste0(index, ".tif"))
  pbapply::pblapply(index, function(i) {
    if (!is.null(x[[i]])) {
      terra::writeRaster(x[[i]], outfiles[i], ...)
    }
  })
  invisible(NULL)
}
