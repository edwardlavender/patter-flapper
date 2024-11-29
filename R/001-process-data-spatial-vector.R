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
library(terra)
library(spatstat)
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

# Define open areas
mpa_open$id     <- seq_len(nrow(mpa_open))
mpa_open$open   <- "open"
mpa_open$col    <- "blue"

# Define closed areas
mpa_closed$id   <- max(mpa_open$id) + 1L
mpa_closed$open <- "closed"
mpa_closed$col  <- scales::alpha("red", 0.5)

# Define mpa SpatVector
mpa <- 
  rbind(mpa_open |> select(id, open, col), 
        mpa_closed |> select(id, open, col)) |> 
  mutate(id = factor(id), 
         open = factor(open)) |>
  st_transform(crs = crs) |> 
  terra::vect()

# Update MPA ids
# 1: Sound of Mull
# 2: South Lismore
# 3: South-East Mull
# 4: Ardnackaig (South Crinan)
# 5: Inverlussa (Jura)
# 6: Luing 
# 7: Crinan
# 8: Scarba 
# 9: Firth of Lorn (closed)
mpa$id <- factor(mpa$id, labels = c("Sound of Mull", "South Lismore", "South-East Mull", 
                                    "Ardnackaig", "Inverlussa", "Luing", 
                                    "Crinan", "Scarba", "Firth of Lorn"))

# Visual checks
data.frame(mpa)
terra::plot(mpa, y = "id", col = 1:9L)
terra::plot(mpa, col = mpa$col)
terra::plot(mpa, y = "open", col = c(mpa_open$col[1], mpa_closed$col[1]))
terra::plot(howe)
terra::lines(mpa, col = scales::alpha(mpa$col, 0.5), lwd = 2)


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
png(here_fig("study-site-quick.png"), 
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

#### Define sea
sea <- st_invert(st_as_sf(coast)) |> terra::vect()
# plot(st_geometry(sea), col = "blue")

#### Define observation window (~4 s)
tic()
win <- as.owin(st_as_sf(sea))
toc()

#### Trial smoother windows (~44 mins)
# dmin = 0      : UD estimation: 161 s
# dmin = 100    : Simplification: 462 s, UD estimation: 19 s
# dmin = 500    : Simplification: 452 s, UD estimation: 12 s
# dmin = 1000   : Simplification: 452 s, UD estimation: 10 s
# dmin = 2000   : Simplification: 442 s, UD estimation: 10 s
# dmin = 5000   : Simplification: 459 s, UD estimation: 10 s
wins <- cl_lapply(c(0, 100, 500, 1000, 2000, 5000), function(dmin) {
  simplify.owin.trial(sea = sea, win = win, dmin = dmin)
})
qs::qsave(wins, here_data("spatial", "wins.qs"))

# > Wiith dmin <= 1000, we have accurate representation of coastline.
# > With dmin = 2000, we start to see the effects of rounding 
# > (without further speedup in UD estimation). 
# > With dmin = 5000, we have too much smoothing of coastline. 
# > (and no further speedup in UD estimation). 

# Use win smoothed by dmin = 1000
win <- wins[[4]] 

#### Visualise observation window 
# (<164 s)
# < 1 s for dmin = 1000
if (FALSE) {
  tic()
  plot(win, col = "blue")
  toc()
}


###########################
###########################
#### Write datasets to file (~1 s)

tic()
terra::writeVector(vect(bb), here_data_raw("boundaries", "bb.shp"))
qsaveext(bb, here_data("spatial", "bb.qs"))
qsavevect(mpa, here_data("spatial", "mpa.qs"))
qsavevect(coast, here_data("spatial", "coast.qs"))
qs::qsave(win, here_data("spatial", "win.qs"))
toc()


#### End of code. 
###########################
###########################