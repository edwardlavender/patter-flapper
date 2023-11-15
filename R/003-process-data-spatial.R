###########################
###########################
#### process-data-spatial.R

#### Aims
# 1) Define study site & associated spatial datasets

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
library(sf)
library(tictoc)

#### Load data
digi       <- terra::rast(here_data_raw("bathymetry", "digmap_bathym_merged_1arcsec_res", "digimap_bathy_merge_1arcsec_res.tif"))
bathy      <- terra::rast(here_data_raw("bathymetry", "Full Data - Original", "EXTRACTED_DEPTH1.tif"))
mpa_open   <- st_read(here_data_raw("mpa", "management_areas", "MPA_Open_area.shp"))
mpa_closed <- st_read(here_data_raw("mpa", "management_areas", "LSSOJ_MPA_CLosed_areat.shp"))
coast      <- st_read(here_data_raw("coast", "westminster_const_region.shp"))


###########################
###########################
#### Process datasets

#### Process coast
coast <- st_transform(coast, terra::crs(bathy))

#### Process bathy (~10 s)
tic()
bathy <- abs(bathy)
names(bathy) <- terra::varnames(bathy) <- "depth"
toc()

#### Collate MPA polygons
# Define open/closed areas
mpa_open$open   <- TRUE
mpa_closed$open <- FALSE
# Join polygons
mpa <- rbind(mpa_open |> select(id, open), 
             mpa_closed |> select(id, open))
# Convert to UTM
mpa <-
  mpa |> 
  st_transform(crs = terra::crs(bathy))

#### Define MPA buffer
# We define an MPA buffer to visualise distances away from the MPA
# We use this as a guide to define a rectangular study area
# * Choose a distance value
# * Visualise the buffer in relation to available bathymetry data
tic()
mpa_buf <- 
  mpa |> 
  st_union() |>
  st_buffer(dist = 30 * 1e3, nQuadSegs = 1e2)
toc()

#### Define study site boundaries based on bathy & mpa_buf
# Define boundary box
terra::ext(bathy)
st_bbox(mpa_buf)
xmin <- 640000
xmax <- 720000
ymin <- 6179750
ymax <- 6318934
# Crop bathy (~10 s)
bathy   <- terra::crop(bathy, c(xmin, xmax, ymin, ymax)) # ~ 10
# Crop MPA (~ 1 s)
mpa_buf <- sf::st_crop(mpa_buf, c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))

#### Merge bathymetry files
# Here, we fill in NAs (e.g., around Lismore) on `bathy` using `digi`
# This is necessary to handle receivers around Lismore
# This introduces some errors:
# * Some NAs will be inappropriately replaced (e.g., near to coastline) 
# * The uncertainty in the bathymetry values varies over space
# * But this enables us to include receivers beyond bathy and, e.g., 
# * ... look at putative movements north of Lismore
tic()
# Mask values on land
mask <- (digi >= 0) + 0
digi <- terra::mask(digi, mask, maskvalues = 1)
# terra::plot(digi)
# terra::plot(is.na(digi))
# - Use absolute values
digi <- abs(digi)
# terra::plot(digi)
# - Align onto bathy geometry
digi <- terra::project(digi, terra::crs(bathy))
digi <- terra::resample(digi, bathy, method = "near")
# Update NAs on bathy with values from digi
mask <- is.na(bathy)
bathy[mask] <- digi[mask]
terra::plot(bathy)
toc()


###########################
###########################
#### Visualise study site 

#### Define graphical parameters
# Define MPA colours for open/closed regions
mpa$col_index <- factor(mpa$open, levels = c("TRUE", "FALSE"))
col_open <- scales::alpha("skyblue", 0.5)
col_closed <- scales::alpha("red", 0.5)
col_mpa <- c(col_open, col_closed)

#### Create plot (~ 1 s)
tic()
png(here_fig("study-site.png"), 
    height = 10, width = 10, units = "in", res = 600)
terra::plot(bathy)
terra::plot(st_geometry(mpa), 
            col = col_mpa[mpa$col_index], 
            lwd = 2, 
            add = TRUE)
if (FALSE) {
  terra::plot(mpa_buf, 
              border = "dimgrey", 
              add = TRUE)
}
terra::lines(coast)
dev.off()
toc()

# > We will create a publication quality map of the study area in QGIS


###########################
###########################
#### Write datasets (~ 15 s)

tic()
saveRDS(mpa, here_data("spatial", "mpa.rds"))
saveRDS(coast, here_data("spatial", "coast.rds"))
terra::writeRaster(bathy, here_data("spatial", "bathy.tif"), overwrite = TRUE)
toc()


#### End of code. 
###########################
###########################