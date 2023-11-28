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
howe <- bathy
mask <- is.na(bathy)
bathy[mask] <- digi[mask]
terra::plot(bathy)
toc()

#### Process coast
coast <- 
  coast |> 
  st_transform(crs = terra::crs(bathy)) |> 
  st_crop(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))


###########################
###########################
#### Analyse bathymetry layers

# Examine the extent of howe on digi
terra::plot(digi)
terra::plot(howe, col = "dimgrey", add = TRUE)
# In regions with data, estimate additional error induced by digi
# * We need to use this information in the depth error model
error <- howe - digi
error <- terra::mask(error, bathy)
# Most errors are small
# * 95 % percentiles: -4.39 m to 7.19 m
terra::hist(error, maxcell = 1e9) 
terra::global(error, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
terra::global(error, "mean", na.rm = TRUE) # 0.5 m (small)
# There are some places where Howe is shallower than Digi
terra::plot(error < -50)
terra::global(error, "min", na.rm = TRUE) 
# There are some areas (coastline) where Howe is MUCH deeper than digi
# * This makes sense given lower resolution of Digi. 
terra::global(error, "max", na.rm = TRUE)  # 349 m (huge!)
# Examine the error with depth 
# * Most errors are 'small'
# * Small errors occur across the full range of depths
# * The largest errors are only found in deep water (this makes sense)
n <- 1e6
mf <- terra::spatSample(error, size = n, 
                        method = "random", replace = FALSE, na.rm = TRUE,
                        values = TRUE, cell = TRUE, as.df = TRUE)
mf$bathy <- terra::extract(howe, mf$cell)[, 1]
colnames(mf) <- c("cell", "error", "howe")
plot(mf$error, mf$howe, pch = ".")
# Examine errors in the region around Lismore
# * This is the most important region where we lack Howe data
# locator()
south_lismore <- terra::ext(705117.5, 719807.8, 
                            6258559, 6271233)
if (FALSE) {
  # Xoom in further:
  south_lismore <- terra::ext(707696.4, 710613.7, 
                              6266425, 6269443)
}
howe_lismore <- terra::crop(howe, south_lismore)
digi_lismore <- terra::crop(digi, south_lismore)
# There is no noticeable jump at the boundary
# * But this map is relatively large in scale
terra::plot(digi_lismore, range = c(0, 250))
terra::plot(howe_lismore, range = c(0, 250), add = TRUE)
terra::sbar(1e3)
# Depths are +/- ~40 m in this region
# * I.e., if an individual moves into this region, 
# * the influence of the depth observations should become weaker
# * That being said, the values are still mostly similar & 
# * ... there is not a strong spatial pattern
error_lismore <- howe_lismore - digi_lismore
terra::global(error_lismore, quantile, na.rm = TRUE)
terra::plot(error_lismore)
terra::plot(error_lismore > 10)
terra::plot(error_lismore < -10)


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
terra::writeRaster(howe, here_data("spatial", "howe.tif"), overwrite = TRUE)
terra::writeRaster(bathy, here_data("spatial", "bathy.tif"), overwrite = TRUE)
toc()


#### End of code. 
###########################
###########################