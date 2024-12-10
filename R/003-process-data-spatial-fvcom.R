###########################
###########################
#### process-data-spatial-fvcom.tbx.R

#### Aims
# (1) Summarise relevant environmental fields in the study area, 
# ... * current speeds, which inform the movement model 
# ... * tidal ranges, which inform the observation model 

#### Prerequisites
# (1) WeStCOMS files (on Lacie_Share)
# (2) Copy WeStCOMS meshes into data/spatial/westcoms (from rzss-flapper)
# Local directory: /Users/lavended/Documents/work/projects/flapper-skate/health/rzss-flapper/rzss-flapper/data/spatial/mesh
# Code: https://github.com/edwardlavender/rzss-flapper/blob/master/R/003_process_spatial.R
# (3) Example code for fvcom.tbx extraction:
# https://github.com/edwardlavender/westcoms_validation/blob/master/R/analyse_study_site.R


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
# try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(dv)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(fvcom.tbx)
library(patter)
library(raster)
library(sp)
library(tictoc)

#### Load data
bathy  <- raster(here_data("spatial", "bathy.tif"))
mesh_around_elements   <- readRDS(here_data("spatial", "westcoms", "mesh_around_elements.rds"))
mesh_around_nodes      <- readRDS(here_data("spatial", "westcoms", "mesh_around_nodes.rds"))


###########################
###########################
#### Preparation

#### Define FVCOM file path
wc_con <- "/Volumes/Lacie_Share/Dima/FVCOM_variable_outputs"
if (!dir.exists(wc_con)) {
  warning("/Volumes/Lacie_Share/ is not connected!")
}

#### Define variable
# variable        <- "current_speed"
variable        <- "tidal_elevation"
wc_con_variable <- file.path(wc_con, variable)
outfile         <- here_data("spatial", "westcoms", paste0(variable, ".qs"))

#### Define mesh
if (variable == "current_speed") {
  mesh <- mesh_around_elements
} 
if (variable == "tidal_elevation") {
  mesh <- mesh_around_nodes
}
# Crop mesh to study area (~26 s)
tic()
crs(mesh) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
mesh      <- spTransform(mesh, crs(bathy))
mesh      <- crop(mesh, extent(bathy))
table(diff(as.numeric(mesh$ID)))
toc()

# (optional) Visualise mesh 
if (FALSE) {
  tic()
  raster::plot(bathy)
  raster::lines(mesh)
  toc()
}

#### Define extract data
# Define parameters
# * We extract data over a one year period
# * Current speeds are extracted for layer 10
# * Data assembly & extraction is implemented iteratively below
dates       <- seq(as.Date("2016-03-01"), as.Date("2017-02-28"), "days")
hours       <- 0:23
layers      <- 10
mesh_IDs    <- as.integer(as.character(mesh$ID))
if (variable == "current_speed") {
  read_fvcom <- readRDS
  ext        <- ".rds"
}
if (variable == "tidal_elevation") {
  read_fvcom <- function(con) R.matlab::readMat(con)$data
  ext        <- ".mat"
}


###########################
###########################
#### Implement extraction

#### Extract model outputs
# Current speed   : ~50 mins (6 cl, 1 cl may be faster!)
# Tidal elevation : ~2 mins (6 cl)
if (!file.exists(outfile)) {
  
  tic()
  wc_by_date <- cl_lapply(dates, .cl = 6L, .fun = function(date) {
    
    # date <- dates[1]
    
    # Define input
    input <- 
      CJ(date_name = date_name(date), 
         hour = hours, 
         layer = layers, 
         mesh_ID = mesh_IDs) |> 
      dplyr::arrange(date_name, 
                     hour, 
                     mesh_ID) |> 
      as.data.frame()
    if (variable == "tidal_elevation") {
      input$layer <- NULL
    }
    
    # Extract information 
    output <- fvcom.tbx::extract(dat = input, # [1:10, ], 
                                 dir2load = wc_con_variable, 
                                 read_fvcom = read_fvcom,
                                 extension = ext, 
                                 verbose = FALSE)
    
    # Return summary
    data.table(n = nrow(output), 
               min = min(output$wc, na.rm = TRUE),
               mean = mean(output$wc, na.rm = TRUE),
               max = max(output$wc, na.rm = TRUE),
               sum = sum(output$wc, na.rm = TRUE), 
               na = length(which(is.na(output$wc))))
    
  })
  toc()
  
  wc <- rbindlist(wc_by_date)
  qs::qsave(wc, outfile)
  
}


###########################
###########################
#### Analysis

stopifnot(all(wc$na == 0))

if (variable == "current_speed") {
 
  # Average, daily-mean current speed
  mean(wc$mean)
  # 0.08595888
  
  # Average current speed over the whole year
  sum(wc$sum) / sum(wc$n)
  # 0.08595888
  # > This is approximately 0.02 DLs/s
  (sum(wc$sum) / sum(wc$n)) / 1.75
  # 0.02148473
  
  # Min/max current speed
  # > These results are as expected (e.g., given high speeds in Gulf of Corryvreckan)
  min(wc$min)
  max(wc$max)
  # 0
  # 3.840138
   
}

if (variable == "tidal_elevation") {
  
  # Min/max tidal elevation 
  min(wc$min)
  max(wc$max)
  # -3.243462
  # 2.659841
  
  # Tidal ranges
  # > Spring tidal ranges are ~1 m
  # > Neap tidal ranges are up to 6 m
  # > These results align with reports
  # > Tidal surges (and waves in shallow water) may elevate effective depths further
  # > An additional 2 m for waves/tidal surges is probably sufficient for the vast majority of the time
  # https://www.gov.scot/publications/scotlands-marine-atlas-information-national-marine-plan/pages/10/# 
  utils.add::basic_stats(wc$max - wc$min)
  # min mean median  max   sd  IQR  MAD
  # 1 1.27 3.38   3.31 5.88 1.07 1.64 1.22
}


#### End of code.
###########################
###########################