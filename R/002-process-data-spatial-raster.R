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
dv::src()

#### Load data
bb    <- qreadext(here_data("spatial", "bb.qs"))
coast <- qreadvect(here_data("spatial", "coast.qs"))
howe  <- terra::rast(here_data_raw("bathymetry", "firth-of-lorn", "EXTRACTED_DEPTH1.tif"))

#### Julia connection
tic()
julia_connect()


###########################
###########################
#### Process Howe dataset (base)

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
terra::plot(howe)
coast |> 
  terra::simplifyGeom(tolerance = 500) |> 
  terra::plot(add = TRUE, border = "dimgrey")


###########################
###########################
#### Expand Howe dataset

# TO DO
bathy <- howe
terra:::readAll(bathy)


###########################
###########################
#### Handle spikes

# Method
# * Iterate over cells (in Julia)
# * Compare the depth of the cell to the median depth for the surrounding eight cells
# * If greater than 50 m below the median, set to the median value
# * TO DO
# * See: https://www.gebco.net/about_us/faq/ 

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

# Identify spikes
spikes <- bathy - smooth
# terra::plot(spikes > 5)


###########################
###########################
#### Validate bathymetry 

# (optional) Compare datasets in overlapping regions
# TO DO
# See old code for outline


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
terra::ncell(bathy)

# Depth range (~2 s)
if (FALSE) {
  tic()
  terra::global(bathy, "range", na.rm = TRUE)
  toc()
  # min     max
  # map_value   0 349.999
}

# Bathymetric uncertainty ranges from 1-5 m 
# > We can envelope this within a single parameter for speed (e.g, 5 m)
ebathy(0)
ebathy(350)

# The total depth of the individual below the seabed may be ~11 m 
# > We will increase this a bit to account for additional uncertainties
# (e.g., multiple datasets)
ebathy(350) + etag + etide


###########################
###########################
#### UD estimation rasters

#### UD grid
# UDs are computed using 500 pixels in the x and y direction
ud_grid <- terra::rast(bb, nrow = 500, ncol = 500, crs = terra::crs(bathy))
ud_grid <- terra::resample(bathy, ud_grid)
ud_grid <- ud_grid > 0
ud_grid <- terra::classify(ud_grid, cbind(0, NA))
ud_grid
terra::plot(!is.na(ud_grid))

#### RSP grid 
# Define 'water' grid
# * See guidance in patter-eval
ud_grid_ll <- terra::project(ud_grid, "EPSG:4326")
names(ud_grid_ll) <- "layer"
terra::plot(ud_grid_ll)
# Define transition layer (~1 s)
tm <- actel::transitionLayer(ud_grid_ll, directions = 8L)
terra::plot(raster::raster(tm))
# Define (enlarged) dynBBMM base.raster
# * Define regional raster from how_full
# * Resample ud_grid onto howe_full for merging
# * Merge the two datasets (i.e., the full dataset with the smaller, refined one)
# * Crop an enlarged dataset (e.g., expanded by 20 km)
tic()
region <- terra::rast(terra::ext(howe_full), nrow = 500, ncol = 500, crs = terra::crs(bathy))
region <- terra::resample(howe_full, region)
bbrast <- terra::resample(ud_grid, region)
bbrast <- terra::merge(bbrast, region, na.rm = FALSE)
bbrast <- terra::crop(bbrast, bb + 20e3L)
bbrast <- bbrast > 0
bbrast <- terra::classify(bbrast, cbind(0, NA))
bbrast_ll <- terra::project(bbrast, "EPSG:4326")
terra::plot(bbrast_ll)
toc()
# Visual check
pp <- par(mfrow = c(1, 3))
terra::plot(howe_full)
terra::plot(ud_grid_ll)
terra::plot(bbrast_ll)
par(pp)


###########################
###########################
#### Write datasets

# Size of bathymetry SpatRaster: 4.5 GB
# * This comprises n doubles
# * Each double is 8 bytes
terra::ncell(bathy)
8 * terra::ncell(bathy) / 1e9L

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
terra::writeRaster(bbrast_ll,
                   here_data("spatial", "bbrast_ll.tif"), 
                   overwrite = TRUE)
qs::qsave(tm, here_data("spatial", "tm.qs"))


toc()


#### End of code. 
###########################
###########################