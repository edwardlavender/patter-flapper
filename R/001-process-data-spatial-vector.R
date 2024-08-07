###########################
###########################
#### process-data-spatial-vector.R

#### Aims
# 1) Define study site & associated spatial (vector) datasets
# * Study site boundaries
# * MPA
# * Coastline

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
howe       <- terra::rast(here_data_raw("bathymetry", "firth-of-lorn", "EXTRACTED_DEPTH1.tif"))
mpa_open   <- st_read(here_data_raw("mpa", "management_areas", "MPA_Open_area.shp"))
mpa_closed <- st_read(here_data_raw("mpa", "management_areas", "LSSOJ_MPA_CLosed_areat.shp"))
coast      <- st_read(here_data_raw("coast", "westminster_const_region.shp"))



###########################
###########################
#### Process coast

# Define CRS
crs <- terra::crs(howe)

# Process coast (~11 s)
tic()
coast <- 
  coast |> 
  st_geometry() |>
  st_transform(crs = crs) |>  
  terra::vect() |> 
  # For downstream speed, crop to Scotland
  terra::crop(terra::ext(c(559885.6, 1069788.1, 6058377, 6437467))) |> 
  terra::aggregate()
toc()
# Plot coast to guess/check coordinates for copping
if (FALSE) {
  coast |> 
    terra::simplifyGeom(tolerance = 5000) |> 
    terra::plot()
  # axis(side = 1)
}


###########################
###########################
#### Process MPA

mpa_open$open   <- TRUE
mpa_open$col    <- scales::alpha("skyblue", 0.5)
mpa_closed$open <- FALSE
mpa_closed$col  <- scales::alpha("red", 0.5)

mpa <- 
  rbind(mpa_open |> select(id, open, col), 
        mpa_closed |> select(id, open, col)) |> 
  st_transform(crs = crs) |> 
  terra::vect()


###########################
###########################
#### Define study area

#### Map approx area
terra::plot(howe)
coast |> 
  terra::simplifyGeom(tolerance = 500) |> 
  terra::plot(add = TRUE, border = "dimgrey")
mpa |> 
  terra::plot(add = TRUE, border = "royalblue", lwd = 2)
mpa |> 
  terra::aggregate() |> 
  terra::buffer(width = 30 * 1e3, quadsegs = 1e2) |> 
  terra::plot(add = TRUE, border = "royalblue", lwd = 2, lty = 2)

#### Define study limits
# locator()
xmin <- 640000
xmax <- 740000 # 720000
ymin <- 6179750
ymax <- 6318934
bb <- terra::ext(c(xmin, xmax, ymin, ymax))
# Expand by 20 km around boundary box (for plotting)
coast <- terra::crop(coast, bb + 100 * 1e3)
# Calculate the size (GB) of a 5 x 5 m bathymetry grid across this region:
ncell <- ((xmax - xmin) / 5) * ((ymax - ymin) / 5)
( gb  <- 8 * ncell / 1e9L )
if (gb > 6) {
  stop("Bathymetry grid > 6 GB in size.")
}

#### Visualise study site (~7 s)
# > We will create a publication quality map of the study area in QGIS
tic()
png(here_fig("study-site.png"), 
    height = 10, width = 10, units = "in", res = 600)
terra::plot(howe)
terra::plot(mpa, 
            col = mpa$open, 
            lwd = 2, 
            add = TRUE)
terra::lines(coast)
dev.off()
toc()


###########################
###########################
#### Prepare {spatstat} inputs

#### Define observation window (~4 s)
tic()
win <- as.owin.sf(st_as_sf(coast), .invert = TRUE)
toc()

#### Visualise observation window (~164 s)
if (FALSE) {
  tic()
  plot(win, col = "blue")
  toc()
}


###########################
###########################
#### Write datasets to file (~1 s)

tic()
qsaveext(bb, here_data("spatial", "bb.qs"))
qsavevect(mpa, here_data("spatial", "mpa.qs"))
qsavevect(coast, here_data("spatial", "coast.qs"))
qs::qsave(win, here_data("spatial", "win.qs"))
toc()


#### End of code. 
###########################
###########################