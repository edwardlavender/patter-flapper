###########################
###########################
#### process-data-spatial-raster.R

#### Aims
# 1) Define bathymetry datasets

#### Prerequisites
# 1) Obtain raw data


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(sf)
library(ggplot2)
dv::src()

#### Load data
bb          <- qreadext(here_data("spatial", "bb.qs"))
coast       <- qreadvect(here_data("spatial", "coast.qs"))
howe        <- terra::rast(here_data_raw("bathymetry", "firth-of-lorn", "EXTRACTED_DEPTH1.tif"))
scotland    <- terra::rast(here_data_raw("bathymetry", "digimap", "scotland", "bathy_scotland_6_arc_sec.tif"))
westminster <- terra::vect(here_data_raw("coast", "westminster_const_region.shp"))

#### Julia connection
tic()
julia_connect()


###########################
###########################
#### Process Howe dataset (base)

if (FALSE) {
  
  # Process base dataset
  howe_full <- run(here_data("spatial", "howe-full.tif"), 
                   overwrite = FALSE, 
                   expr = {
                     # (~33 s)
                     howe <- abs(howe)
                     names(howe) <- terra::varnames(howe) <- "map_value"
                     howe
                   }, 
                   read = terra::rast, 
                   write = terra::writeRaster)
  
  # Visualise full dataset & gaps:
  # * Loch Linnhe
  # * South of Jura
  # * West of Islay/Coll & Tiree
  # * South and Eastern Skye
  # * West and Northern Skye
  terra::plot(howe_full)
  coast |> 
    terra::simplifyGeom(tolerance = 500) |> 
    terra::plot(add = TRUE, border = "dimgrey")
  
  # Crop howe 
  howe <- run(here_data("spatial", "howe.tif"), 
              overwrite = FALSE, 
              expr = {
                # ~38 s
                terra::crop(howe_full, bb)
              }, 
              read = terra::rast, 
              write = terra::writeRaster)
  
  # Visualise cropped dataset & gaps:
  # * Loch Linnhe
  # * North-western Coll
  # * Southern boundary (South Jura)
  plot_howe <- function() {
    terra::plot(howe)
    coast |> 
      terra::simplifyGeom(tolerance = 500) |> 
      terra::plot(add = TRUE, border = "dimgrey")
    nothing()
  }
  plot_howe()
  
}


###########################
###########################
#### Expand Howe dataset

if (FALSE) {
  
  bathy <- howe
  readAll(bathy)
  
  #### Loch Linnhe (~36 s)
  files <- list.files(here_data_raw("bathymetry", "loch-linnhe-and-etive", "datasets"), 
                      pattern = "Loch Linnhe Blk.*\\.bag$", 
                      full.names = TRUE, recursive = TRUE)
  files_linnhe <- files
  terra::rast(files[1])
  linnhe <- get_bathy(files = files, 
                      bathy = bathy, 
                      coast = coast, 
                      outfile = here_data_raw("bathymetry", "temporary", "tiles", "loch-linnhe.tif"))
  
  #### Loch Etive (~11 s)
  files <- list.files(here_data_raw("bathymetry", "loch-linnhe-and-etive", "datasets"), 
                      pattern = "Loch Etive.*\\.bag$", 
                      full.names = TRUE, recursive = TRUE)
  files_etive <- files
  terra::rast(files[1])
  etive <- get_bathy(files = files,
                     bathy = bathy, 
                     coast = coast, 
                     outfile = here_data_raw("bathymetry", "temporary", "tiles", "loch-etive.tif"))
  terra::plot(linnhe, add = TRUE, col = scales::alpha("purple", 0.5))
  
  #### Loch Creran (~6 s)
  files <- list.files(here_data_raw("bathymetry", "loch-creran", "datasets"), 
                      pattern = "Loch Creran.*\\.bag$", 
                      full.names = TRUE, recursive = TRUE)
  files_creran <- files
  terra::rast(files[1])
  creran <- get_bathy(files = files,
                      bathy = bathy, 
                      coast = coast, 
                      outfile = here_data_raw("bathymetry", "temporary", "tiles", "loch-creran.tif"))
  terra::plot(linnhe, add = TRUE, col = scales::alpha("purple", 0.5))
  
  #### Coll and Tiree (~82 s)
  files <- list.files(here_data_raw("bathymetry", "coll-and-tiree", "datasets"), 
                      pattern = "\\.bag$", 
                      full.names = TRUE, recursive = TRUE)
  files_coll <- files
  terra::rast(files[1])
  coll <- get_bathy(files = files,
                    bathy = bathy, 
                    coast = coast, 
                    outfile = here_data_raw("bathymetry", "temporary", "tiles", "col.tif"))
  
  #### Islay (~34 mins)
  files <- list.files(here_data_raw("bathymetry", "islay", "datasets"), 
                      pattern = "\\.bag$", 
                      full.names = TRUE, recursive = TRUE)
  files_islay <- files
  terra::rast(files[1])
  islay <- get_bathy(files = files,
                     bathy = bathy, 
                     coast = coast, 
                     outfile = here_data_raw("bathymetry", "temporary", "tiles", "islay.tif"))
  
  #### Lochgilphead (~89 s)
  files <- list.files(here_data_raw("bathymetry", "lochgilphead", "datasets"), 
                      pattern = "\\.bag$", 
                      full.names = TRUE, recursive = TRUE)
  files_lochgilphead <- files
  terra::rast(files[1])
  lochgilphead <- get_bathy(files = files,
                            bathy = bathy, 
                            coast = coast, 
                            outfile = here_data_raw("bathymetry", "temporary", "tiles", "lochgilphead.tif"))
  
  #### Merge datasets (~3 mins?)
  tic()
  # Expand datasets
  bathys <- list(linnhe, etive, creran, coll, islay, lochgilphead)
  bathys <- lapply(bathys, function(r) {
    terra::resample(r, bathy, threads = TRUE)
  })
  # Merge datasets (~2 mins)
  bathys <- append(list(bathy), bathys)
  bathy <- do.call(terra::merge, bathys)
  # Save dataset
  terra::writeRaster(bathy, 
                     here_data_raw("bathymetry", "temporary", "merged-dataset.tif"), 
                     overwrite = TRUE)
  toc()
  
  #### Record datasets
  dlinnhe       <- files_linnhe |> sort_bathysets()
  detive        <- files_etive |> sort_bathysets()
  dcreran       <- files_creran |> sort_bathysets()
  dcoll         <- files_coll |> sort_bathysets()
  dislay        <- files_islay |> sort_bathysets()
  dlochgilphead <- files_lochgilphead |> sort_bathysets()
  bathysets <- data.table(
    Region = c("(Base) Firth of Lorn", 
               fill(dcoll, "Coll and Tiree"),
               fill(dislay, "Islay"), 
               fill(dcreran, "Loch Creran"), 
               fill(detive, "Loch Etive"), 
               fill(dlinnhe, "Loch Linnhe"), 
               fill(dlochgilphead, "Lochgilphead")),
    Dataset = c("Howe et al. (2014)", dcoll, dislay, dcreran, detive, dlinnhe, dlochgilphead))
  prettyGraphics::tidy_write(bathysets, here_fig("bathysets.txt"))
  # rstudioapi::navigateToFile(here_fig("bathysets.txt"))
  
} 

bathy <- terra::rast(here_data_raw("bathymetry", "temporary", "merged-dataset.tif"))


###########################
###########################
#### Handle inshore areas

# (optional) TO DO

# There are inshore areas where bathymetry data are lacking 
# E.g., around Lismore where data for receiver 2 are lacking 
# We could fill these data with zeros
# This gives the algorithms more 'shallow-water' space to explore
# But coastline data are also uncertain
# This is not currently implemented

# Method:
# Define grid with zeros
# Mask areas in the sea (coast + 100 m buffeer)
# Mask areas on land (coast)
# Keep zeros in a 100 m band around coastline
# Merge with bathy, replacing NAs in this region with zeros, if required


###########################
###########################
#### Handle spikes

# Method
# * Iterate over cells (in Julia)
# * Compare the depth of the cell to the median depth for the surrounding eight cells
# * If greater than 50 m below the median, set to the median value
# * For further information about spikes, see https://www.gebco.net/about_us/faq/ 

if (FALSE) {
  
  # Compute the median value for each cell based on surrounding cells (~46 s)
  tic()
  focals <- terra::focal(bathy, w = 3, 
                         fun = "median", na.rm = TRUE,
                         na.policy = "omit")
  toc()
  
  # Export data to Julia for processing (~15 s)
  tic()
  julia_assign_SpatRaster("bathy", bathy)
  julia_assign_SpatRaster("focals", focals)
  toc()
  
  # Process spikes in Julia (~20 s)
  tic()
  julia_code(
    '
  using Base.Threads: @threads
  @time for i in 1:size(bathy)[1]
    println(i)
    @threads for j in 1:size(bathy)[2]
      if bathy[i, j] > (focals[i, j] + 50)
        bathy[i, j] = focals[i, j]
      end
    end
  end
  GeoArrays.write("data/spatial/spikes.tif", bathy)
  '
  )
  toc()
  
}

# Load smoothed bathymetry into R
smooth <- terra::rast("data/spatial/spikes.tif")
names(smooth) <- "map_value"

# Identify spikes
spikes <- bathy - smooth
# terra::plot(spikes > 5)

# Use smoothed bathymetry 
bathy <- smooth


###########################
###########################
#### Validate bathymetry 

# (optional) Compare datasets in overlapping regions
# (optional) TO DO
# See old code for outline


###########################
###########################
#### Aggregation

#### Motivation
# With a high-resolution bathymetry layer convergence is challenging
# Very large numbers of particles are required to achieve convergence, even in simulations, which is expensive
# An aggregated bathymetry layer may facilitate convergence and reduce computation time
# Here, we explore aggregation options and the induced error

#### Record bathy 5 m
bathy_5m <- terra::deepcopy(bathy)
rnow     <- terra::res(bathy_5m)
terra::writeRaster(bathy_5m, here_data("spatial", "bathy-5m.tif"), overwrite = TRUE)

#### Aggregate bathymetry 
# Define desired resolution
ud_grid <- terra::rast(bb, nrow = 500, ncol = 500, crs = terra::crs(bathy))
rnew    <- terra::res(ud_grid)
rnew    <- c(rnew[2], rnew[1])
# We use a ~500 x 500 pixel raster for UD estimation
bathy_agg <- terra::aggregate(bathy, fact = rnew / rnow, fun = "mean", na.rm = TRUE)
# Visualise aggregation
# > There are very few spikes around the coastline (good)
terra::plot(bathy_agg > 100)

#### Quantify bathymetric variation within aggregated cells
# Quantify range within aggregated cells
bathy_agg_max <-  terra::aggregate(bathy, fact = rnew / rnow, fun = "max", na.rm = TRUE)
bathy_agg_min <- terra::aggregate(bathy, fact = rnew / rnow, fun = "min", na.rm = TRUE)
bathy_agg_rng <- bathy_agg_max - bathy_agg_min
terra::hist(bathy_agg_rng)
terra::global(bathy_agg_rng, "sd", na.rm = TRUE)
terra::global(bathy_agg_rng, quantile, probs = seq(0, 1, by = 0.01), na.rm = TRUE) 
# Check SD within aggregated cells 
# > The max SD reaches ~152.3439 m but this may be due to spikes
# > We use a more conservative value half this magnitude (75 m)
# > This corresponds to ~99.95 percentile & covers the vast majority of cells
map_agg_sd <- terra::aggregate(bathy, fact = rnew / rnow, fun = "sd", na.rm = TRUE)
terra::global(map_agg_sd, "mean", na.rm = TRUE)
terra::global(map_agg_sd, "max", na.rm = TRUE)
terra::global(map_agg_sd, quantile, prob = seq(0.95, 1, by = 0.001), na.rm = TRUE)
terra::global(map_agg_sd, quantile, prob = 0.9995, na.rm = TRUE) 
# Visualise SD
# * This is high along coastlines (steep relief)
terra::plot(map_agg_sd)
# Check maximum depth below depth of aggregated cell (+ 20 m)
# * The differences are greatest near coastline
bathy_agg_max   <- terra::aggregate(bathy, fact = rnew / rnow, fun = "max", na.rm = TRUE)
bathy_agg_delta <- bathy_agg_max - bathy_agg
terra::plot(bathy_agg_delta)
terra::global(bathy_agg_delta, quantile, prob = seq(0.9, 1, by = 0.001), na.rm = TRUE)

#### Aggregate bathy & compute error in terms of differences in depth
# Compute aggregation error for different levels of aggregation (~7 mins)
# c(10, 25, etc.) are potential choices for the resolution of the bathymetry raster
if (FALSE) {
  
  aggerror <- 
    cl_lapply(list(10, 25, 50, 75, 100, 200, 250), function(rnew) {
      
      # Aggregate bathy to selected resolution 
      # ~4s for 50 m resolution
      tic()
      bathy_agg <- terra::aggregate(bathy, fact = rnew / rnow, fun = "max", na.rm = TRUE)
      toc()
      
      # Disaggregate bathy to original resolution 
      # ~16 s for 50 to 5 m 
      tic()
      bathy_disagg <- terra::resample(bathy_agg, bathy_5m, method = "near", threads = TRUE) 
      toc()
      
      # Compute error induced by aggregation
      bathy_delta <- bathy_disagg - bathy_5m
      
      # Visually examine error
      # * Note there are some v. large errors as not all spikes have been processed
      # terra::plot(bathy_delta)
      
      # Compute error statistics (~30 s)
      # terra::global(bathy_delta, "mean", na.rm = TRUE)    # 2.275835
      # terra::global(bathy_delta, median, na.rm = TRUE)    # 0.6148987 
      # terra::global(bathy_delta, "sd", na.rm = TRUE)      # 10.9404
      tic()
      probs  <- seq(0.9, 1, by = 0.01)
      errors <- terra::global(bathy_delta, quantile, probs = probs, na.rm = TRUE)
      toc()
      
      # Return data.table
      data.table(rnew = rnew, probs = probs, error = as.numeric(errors[1, ]))
      
    }) |> rbindlist()
  
  #### Visualise error ~ resolution (by probs quantile)
  # The error increases slightly with at aggregated resolutions
  # This is similar except for 0.98 & 0.99 quantiles where the error grows more quickly
  png(here_fig("bathy-aggregation-error.png"), 
      height = 8, width = 10, units = "in", res = 600)
  aggerror |>
    ggplot(aes(rnew, error)) + 
    geom_line() + 
    facet_wrap(probs) |> 
    print()
  dev.off()
  
  #### Aggregate bathymetry
  aggerror[rnew == 100, ]
  # rnew probs      error
  # <num> <num>      <num>
  #   1:   100  0.90   8.747966
  # 2:   100  0.91   9.473255
  # 3:   100  0.92  10.327118
  # 4:   100  0.93  11.357324
  # 5:   100  0.94  12.629799
  # 6:   100  0.95  14.281601
  # 7:   100  0.96  16.560842
  # 8:   100  0.97  20.054390
  # 9:   100  0.98  26.463852
  # 10:   100  0.99  53.530594
  # 11:   100  1.00 349.974457
  # rnew      <- 100
  # bathy_agg <- terra::aggregate(bathy, fact = rnew / rnow, fun = "max", na.rm = TRUE)
  # terra::plot(bathy_agg)
  
}

#### Use aggregated bathymetry layer
bathy_agg    <- terra::resample(bathy_agg, ud_grid, method = "bilinear", threads = TRUE)
bathy        <- terra::deepcopy(bathy_agg)
names(bathy) <- "map_value"
bathy
# terra::plot(bathy)

#### Map bbox
bathy_bbox <- patter:::map_bbox(bathy)


###########################
###########################
#### Bathymetry statistics 

# Study area size (~38 s)
if (FALSE) {
  tic()
  area <- terra::cellSize(bathy, unit = "km")
  terra::global(area, "sum")
  # 13916.83
  toc()
}

# Number of cells
terra::ncell(bathy_5m) # 556740000
terra::ncell(bathy)    # 253500

# Depth range (~2 s)
if (FALSE) {
  tic()
  terra::global(bathy_5m, "range", na.rm = TRUE)
  toc()
  # min     max
  # map_value   0 349.999
}

# Bathymetric uncertainty for the high-resolution dataset
if (FALSE) {
  
  # Bathymetric uncertainty ranges from 1-5 m 
  # > We can envelope this within a single parameter for speed (e.g, 5 m)
  ebathy(0)
  ebathy(350)
  
  # The total depth of the individual below the seabed may be ~11 m 
  # > We will increase this a bit to account for additional uncertainties
  # (e.g., multiple datasets)
  ebathy(350) + etag + etide
  
  # Uncertainties for Order 1a surveys and other types
  # https://iho.int/uploads/user/pubs/standards/s-44/S-44_Edition_6.1.0.pdf, page 18
  # * Howe et al. (2014) survey was order 1a (confirmed by John Howe)
  # * Order 1b surveys are the same in terms of TVU
  # * Order 2 surveys are less accurate (± 8 m): 
  # a = 1.0 m
  # b = 0.023
  ebathy(0, .a = 1.0, .b = 0.023)
  ebathy(350, .a = 0.5, .b = 0.023)
  
}


###########################
###########################
#### UD estimation rasters

#### UD grid
# UDs are computed using 500 pixels in the x and y direction
ud_grid <- terra::resample(bathy, ud_grid)
ud_grid <- ud_grid > 0
ud_grid <- terra::classify(ud_grid, cbind(0, NA))
ud_grid
terra::plot(!is.na(ud_grid))

#### UD grid (NULL model)
ud_null <- terra::setValues(ud_grid, 1)
ud_null <- terra::mask(ud_null, ud_grid)
ud_null <- spatNormalise(ud_null)

#### RSP grid 
# For RSPs, we use an enlarged bathymetry layer
# 1) Cut bathymetry layer by coastline (~3 s)
# westminster <- terra::simplifyGeom(westminster, tolerance = 100)
tic()
westminster <- terra::project(westminster, terra::crs(scotland))
scotland    <- terra::mask(scotland, westminster, inverse = TRUE)
toc()
# 2) Define 'water'
scotland <- abs(scotland)
scotland <- scotland > 0
scotland <- terra::classify(scotland, cbind(0, NA))
# 3) Visual check 
scotland    <- terra::project(scotland, terra::crs(ud_grid), threads = TRUE)
westminster <- terra::project(westminster, terra::crs(ud_grid))
terra::plot(scotland)
# terra::lines(westminster)
terra::lines(bb, col = "red")
bbrsp <- bb + 200000
terra::lines(bbrsp, col = "red")
# 4) Crop bathymetry layer to sensible region 
scotland <- terra::crop(scotland, bbrsp)
# 5) Resample to match ud_grid 
scotland_blank <- terra::rast(bbrsp, res = terra::res(ud_grid))
scotland       <- terra::resample(scotland, scotland_blank)
# 6) Force ud_grid and the ud_grid for RSPs to match in the overlapping region
ud_grid_tmp <- terra::deepcopy(ud_grid)
terra::origin(ud_grid_tmp) <- terra::origin(scotland)
scotland <- terra::merge(ud_grid_tmp, scotland)
terra::plot(scotland)
# 7) Convert to WGS84 for RSPs 
ud_grid_ll <- terra::project(scotland, "WGS84", threads = TRUE)
names(ud_grid_ll) <- "layer"
terra::plot(ud_grid_ll)
# 8) Define transition layer (~42 s)
tic()
tm <- actel::transitionLayer(ud_grid_ll, directions = 8L)
terra::plot(raster::raster(tm))
toc()
# 9) Rename to bbrast_ll (backwards compatibility)
bbrast_ll <- terra::deepcopy(ud_grid_ll)
terra::plot(bbrast_ll)
# terra::lines(bathy |> terra::project("EPSG:4326", threads = TRUE) |> terra::ext())


###########################
###########################
#### Write datasets

# Size of full bathymetry SpatRaster: 4.5 GB
# * This comprises n doubles
# * Each double is 8 bytes
terra::ncell(bathy_5m)
8 * terra::ncell(bathy_5m) / 1e9L

# Size of aggregated bathymetry SpatRaster: 2 MB
terra::ncell(bathy)
8 * terra::ncell(bathy) / 1e6

# Size of UD rasters: < 1 MB
8 * terra::ncell(ud_grid) / 1e9L
8 * terra::ncell(bbrast_ll) / 1e9L

# Write data
terra::writeRaster(bathy,
                   here_data("spatial", "bathy.tif"), 
                   overwrite = TRUE)
terra::writeRaster(ud_grid,
                   here_data("spatial", "ud-grid.tif"), 
                   overwrite = TRUE)
terra::writeRaster(ud_null,
                   here_data("spatial", "ud-null.tif"), 
                   overwrite = TRUE)
terra::writeRaster(bbrast_ll,
                   here_data("spatial", "bbrast_ll.tif"), 
                   overwrite = TRUE)
qs::qsave(bathy_bbox, here_data("spatial", "bathy-bbox.qs"))
qs::qsave(tm, here_data("spatial", "tm.qs"))

toc()


#### End of code. 
###########################
###########################
