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
library(dv)
library(dplyr)
library(JuliaCall)
library(patter)
library(sf)
library(tictoc)
dv::src()

#### Load data
bb    <- qreadext(here_data("spatial", "bb.qs"))
coast <- qreadvect(here_data("spatial", "coast.qs"))
howe  <- terra::rast(here_data_raw("bathymetry", "firth-of-lorn", "EXTRACTED_DEPTH1.tif"))

#### Julia connection
julia_connect()


###########################
###########################
#### Process Howe dataset (base)

# Process base dataset
howe <- run(here_data("spatial", "howe-full.tif"), 
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
terra::plot(howe)
coast |> 
  terra::simplifyGeom(tolerance = 500) |> 
  terra::plot(add = TRUE, border = "dimgrey")

# Crop howe 
howe <- run(here_data("spatial", "howe.tif"), 
            overwrite = FALSE, 
            expr = {
              # ~38 s
              terra::crop(howe, bb)
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
#### Bathymetric uncertainty 

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
#### Write datasets

# Size of bathymetry SpatRaster: 3.5 GB
# * This comprises n doubles
# * Each double is 8 bytes
terra::ncell(bathy)
8 * terra::ncell(bathy) / 1e9L

# Write data
terra::writeRaster(bathy,
                   here_data("spatial", "bathy.tif"), 
                   overwrite = TRUE)

# UD grid (25 x 25 m)
ud_grid <- terra::aggregate(bathy, fact = 5)
terra::writeRaster(ud_grid,
                   here_data("spatial", "ud-grid.tif"), 
                   overwrite = TRUE)



#### End of code. 
###########################
###########################