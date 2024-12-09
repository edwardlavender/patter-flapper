if (!patter:::os_linux() | Sys.getenv("JULIA_SESSION") == "FALSE") {

get_bathy <- function(files, bathy, coast, outfile = NULL) {
  
  tictoc::tic()
  on.exit(tictoc::toc(), add = TRUE)
  
  # Read rasters, project onto UTM 29N
  region_rasts_raw <- cl_lapply(files, function(x) {
    r <- terra::rast(x)
    if ("elevation" %in% names(r)) {
      r <- r$elevation
    }
    if (terra::nlyr(r) > 1L) {
      stop("Raster contains multiple layers.")
    }
    terra::project(r, terra::crs(bathy), threads = TRUE)
  })
  
  # Check for uncertainty
  # * For Loch Linnhe, this is not available
  pp <- one_page(TRUE, length(files))
  cl_lapply(files, function(x) {
    r <- terra::rast(x)
    # print(names(r))
    if (!("uncertainty" %in% names(r))) {
      warning("uncertainty layer is absent")
    }
    r <- r$uncertainty
    terra::plot(r)
  })
  par(pp)
  
  # Determine extent 
  region_bbs <- do.call(rbind, 
                        lapply(region_rasts_raw, \(x) {
                          e <- terra::ext(x)
                          c(e[1], e[2], e[3], e[4])
                        })
  )
  region_bb <- c(min(region_bbs[, 1]), max(region_bbs[, 2]), 
                 min(region_bbs[, 3]), max(region_bbs[, 4]))
  region_bb <- terra::intersect(terra::ext(bathy), terra::ext(region_bb))
  
  # Define blank raster for region (for resampling -> merger)
  region_blank <- terra::crop(bathy, region_bb)
  
  # Resample rasters for region onto blank map
  # * This is required for terra::merge()
  region_rasts <- lapply(region_rasts_raw, function(r) {
    # r <- terra::crop(r, region_bb)
    r <- abs(r)
    terra::resample(r, region_blank, threads = TRUE)
  })
  
  # Merge resampled rasters 
  if (length(region_rasts) == 1L) {
    region <- region_rasts[[1]]
  } else {
    region <- do.call(terra::merge, region_rasts)
  }
  
  # Visual check
  terra::plot(bathy)
  coast |> 
    terra::simplifyGeom(tolerance = 500) |> 
    terra::plot(add = TRUE, border = "dimgrey")
  terra::plot(region, col = scales::alpha("red", 0.5), add = TRUE)

  # Write to file
  if (!is.null(outfile)) {
    terra::writeRaster(region, 
                       outfile, 
                       overwrite = TRUE)
  }

  # Return map
  region 
  
}

}
