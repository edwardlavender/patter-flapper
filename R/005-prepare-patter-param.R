###########################
###########################
#### prepare-patter-param.R

#### Aims
# 1) Prepare patter parameters (e.g., mobility)

#### Prerequisites
# 1) flapper_appl project


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(dv)
library(data.table)
library(dtplyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(collapse)
library(patter)
library(tictoc)

#### Load data
dv::src()
howe      <- terra::rast(here_data("spatial", "howe.tif"))
bathy     <- terra::rast(here_data("spatial", "bathy.tif"))
moorings  <- readRDS(here_data("mefs", "moorings.rds"))
acoustics <- readRDS(here_data("mefs", "acoustics.rds"))
archival  <- readRDS(here_data("mefs", "archival.rds"))
pars      <- readRDS(here_data("input", "pars.rds"))

###########################
###########################
#### Define parameters

mobility         <- pars$flapper$mobility
gamma            <- pars$flapper$detection_range
calc_depth_error <- pars$flapper$calc_depth_error


###########################
###########################
#### False detections

# There are no false detections flagged by glatos
# This code is commented to avoid {glatos} installation via {renv}
# acoustics |> 
#   as_glatos() |> 
#   glatos::false_detections(tf = 30 * 60)


###########################
###########################
#### Detection time stamps

#### Validate detection timing 
# > The time between sequential detections should never be < 30 s
acoustics |> 
  group_by(individual_id) |> 
  arrange(timestamp) |>
  summarise(min_gap = min(min(serial_difference(timestamp), na.rm = TRUE))) |>
  arrange(min_gap) |>
  filter(min_gap != 0) |> 
  filter(min_gap < 30) |>
  as.data.table()

# There are gaps < 30 s between sequential detections
# This may be due to clock drift
# We aggregate time stamps to 2 minute intervals to mitigate this issue


###########################
###########################
#### Mobility: from acoustics

# We have previously estimated mobility as 500 m per 2 mins 
# Here, we confirm that there are no transitions between receivers exceeding mobility 

# Define acoustic detections with associated coordinates
acoustics <- 
  acoustics |> 
  left_join(moorings |> 
              select("receiver_id", "receiver_easting", "receiver_northing", "receiver_range"), 
            by = "receiver_id")

# Calculate speeds 
acoustics <- 
  acoustics |> 
  group_by(individual_id) |> 
  arrange(timestamp) |> 
  mutate(
    # Calculate minimum travel distance (between nearest detection container edges)
    dist = dist_along_path(cbind(receiver_easting, receiver_northing), 
                           .lonlat = FALSE) - gamma*2, 
    # Translate distances to speeds (metres per minute)
    duration = serial_difference(timestamp, units = "mins"), 
    speed = dist / as.numeric(duration)) |> 
  as.data.table()

# Validate there are no movements exceeding mobility 
table(acoustics$speed > (mobility / 2))


###########################
###########################
#### Mobility: from archival 


###########################
#### Vertical movement

#### Vertical distances
# (ignore breaks in time series)
tic()
archival <- 
  archival |> 
  group_by(individual_id) |>
  mutate(vdist = abs(serial_difference(depth))) |> 
  filter(!is.na(vdist)) |>
  ungroup() |> 
  as.data.table()
archival |>
  ggplot(aes(vdist)) + 
  geom_histogram() +
  # geom_rug() + 
  facet_wrap(~individual_id)
toc()

#### Autocorrelation
# There is strong autocorrelation in depth time series.
# This suggests strong autocorrelation in horizontal distances, 
# (but there is not a strong relationship between vertical & horizontal distances, 
# see below).

# Example time series:
arc <- archival[individual_id == 6, ]
plot(arc$timestamp, arc$depth, type = "l")

# Full ACF
# * Note the periodic patterms
rho <- acf(arc$depth, nrow(arc))
rho$acf[2]
# * ~1000 time steps for a reversal
plot(rho, xlim = c(0, 3000))

# Simulation ACF 
# * ACF for simulated data is qualitatively similar (but distinct)
ar1_ts <- arima.sim(n = nrow(arc), list(ar = rho$acf[2]))
acf(ar1_ts, nrow(arc))
acf(ar1_ts, 3000)

# Autocorrelation as a function of vertical movement
# * Vertical movements are less autocorrelated than depth observations
# * Small vertical distances are more autocorrelated than large vertical distances
# * (on average)
rho  <- acf(arc$vdist, nrow(arc))
mrho <- data.frame(rho = rho$acf, 
                   vdist = arc$vdist)
mrho |>
  ggplot(aes(vdist, rho)) + 
  geom_point()


###########################
#### Vertical and horizontal movement

# Approach:
# Can we improve upon our estimates of mobility/movement parameters by using 
# the relationship between bwtn vertical dists and horiz. dists
# (such that we can use known vertical distances, from archival time series, 
# to learn about horizontal travel distances in the study area, assuming
# skate behave benthically.)

# Method:
# * Sample pairs of points in the study area;
# * Compare vertical and horizontal distances;
# * Expect shorter vertical distances are generally (but not exclusively)
# * ... associated with shorter horizontal (travel) distances.

# Summary results:
# * This is not really what we find
# * It is difficult to see a relationship between vertical and horizontal dists
# * Small vertical distances may be associated with short/long horiz. dists
# * Large vertical distances may be associated with short/long horiz. dists
# * TLDL: we can't use archival time series as a good indication of mobility

#### Simulate transects
tic()
n   <- 1e6L
xy0 <- terra::spatSample(bathy, size = n, xy = TRUE, na.rm = TRUE)
len <- runif(n, terra::res(bathy)[1], pars$patter$mobility)
ang <- rwn(n)
xy1 <- cstep(.xy0 = xy0[, c("x", "y")], 
             .len = len, 
             .ang = ang, 
             .lonlat = FALSE) 
tic()

#### Collate transect depths and distances
comp <- data.table(
  x0 = xy0$x, y0 = xy0$y, z0 = xy0$depth, 
  x1 = xy1[, 1], y1 = xy1[, 2], z1 = terra::extract(bathy, xy1)$depth, 
  hdist = len)
# Drop points on land 
comp <- comp[!is.na(z1), ]
# Calculate vertical distances (m)
comp[, vdist := abs(z0 - z1)]
# (optional) Define bins for plotting
# * We visualise below how the relationship between vdist and hdist changes
# * ... at different vdist values
bw <- 50
comp[, bin := cut(vdist, breaks = seq(min(vdist), max(vdist) + bw, by = bw), 
                  right = FALSE)]
nrow(comp)

#### Visualise the relationships between vertical and horizontal distances
# Scatter plot (~15 s)
png(here_fig("vdist-vs-hdist-scatter.png"), 
    height = 10, width = 10, units = "in", res = 600)
tic()
comp |>
  ggplot(aes(vdist, hdist)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  xlim(range(comp$vdist)) + 
  ylim(range(comp$hdist))
toc()
dev.off()
# Scatter plot, for different distance bands (~15 s)
png(here_fig("vdist-vs-hdist-scatter-facets.png"), 
    height = 10, width = 10, units = "in", res = 600)
tic()
comp |>
  ggplot(aes(vdist, hdist)) + 
  geom_point() + 
  facet_wrap(~bin)
toc()
dev.off()
# 2d histogram (~1 s)
png(here_fig("vdist-vs-hdist-hist.png"), 
    height = 10, width = 10, units = "in", res = 600)
tic()
comp |>
  ggplot(aes(vdist, hdist)) + 
  geom_bin2d(binwidth = c(50, 50))
toc()
dev.off()
# 2d histogram, for different distance bands (~1 s)
# * Use lapply() so each panel has an independent colour scheme
png(here_fig("vdist-vs-hdist-hist-facets.png"), 
    height = 10, width = 10, units = "in", res = 600)
tic()
cowplot::plot_grid(
  plotlist = 
    lapply(split(comp, comp$bin), function(d) {
      d |>
        ggplot(aes(vdist, hdist)) + 
        geom_bin2d(binwidth = c(5, 20)) + 
        xlim(c(0, 350)) + 
        ggtitle(d$bin[1])
    }))
toc()
dev.off()

#### Model
# mod <- lm(hdist ~ 0 + vdist, data = comp)

#### Results
# * At small vertical distances (5 m), horizontal distances range from 0 - 500 m
# * (But smaller horizontal distances are more common)
# * At larger vertical distances, horizontal distances range from 0 - 500 m
# * (But larger vertical distances are more common)
# * But this relationship is highly variable.


###########################
###########################
#### Detection containers

# Validate that all ~simultaneous detections occur at overlapping receivers

#### Define simultaneous detections at multiple receivers
overlapping <- 
  acoustics |>
  # Round time steps 
  mutate(timestamp = lubridate::round_date(timestamp, "2 mins")) |>
  group_by(individual_id, receiver_id, timestamp) |>
  slice(1L) |>
  ungroup() |>
  # List receivers with detections at each time step
  group_by(individual_id, timestamp) |>
  arrange(receiver_id, .by_group = TRUE) |>
  summarise(n = n(), 
            receiver_id = list(unique(receiver_id)), 
            receiver_key = paste(receiver_id, collapse = ", ")) |> 
  ungroup() |>
  # Identify time steps with detections at n > 1 receiver
  select(receiver_id, receiver_key, n) |> 
  filter(n > 1L) |>
  filter(!duplicated(receiver_key)) |> 
  mutate(index = row_number()) |>
  as.data.table() 

#### Validate container overlap (interactively)
# > All simultaneous detections are at overlapping receivers 
# > (for gamma = 750 m)
if (FALSE) {
  lapply(split(overlapping, seq_row(overlapping)), function(d) {
    
    # Define receiver coordinates
    # d <- overlapping[27, ]
    d <- 
      d |> 
      tidyr::unnest(cols = "receiver_id") |> 
      left_join(moorings |> 
                  filter(receiver_id %in% receiver_id), 
                by = "receiver_id") |> 
      as.data.table()
    
    # Define container & plot 
    container <- 
      lapply(seq_row(d), function(i) {
        b <- 
          cbind(d$receiver_easting[i], d$receiver_northing[i]) |> 
          terra::vect() |> 
          terra::buffer(width = gamma)
        if (i == 1L) {
          terra::plot(b, 
                      xlim = c(min(d$receiver_easting - 1e3), max(d$receiver_easting) + 1e3),
                      ylim = c(min(d$receiver_northing - 1e3), max(d$receiver_northing) + 1e3), 
                      main = nrow(d))
        } else {
          terra::lines(b)
        }
        b
      }) |> 
      patter:::spatIntersect()
    terra::lines(container, lwd = 4, lty = 3, col = "red")
    text(d$receiver_easting, d$receiver_northing, d$receiver_id, font = 2)
    
    # Continue 
    readline("Press [Enter] to continue...")
    NULL
  }) |> invisible()
}


###########################
###########################
#### Concurrent acoustic & depth observations 

# At the moment of detection, we know the area within which an individual must be located
# Here, we validate there are possible locations within each _detection_ container
# & tune our gamma & depth error model parameters accordingly 

#### Isolate acoustics data with associated archival observations
acoustics_w_arc <- 
  acoustics |> 
  filter(individual_id %in% archival$individual_id) |> 
  group_by(individual_id) |> 
  filter(timestamp %within% 
           lubridate::interval(
             min(archival$timestamp[archival$individual_id == individual_id[1]]),
             max(archival$timestamp[archival$individual_id == individual_id[1]]))
         ) |> 
  ungroup() |>
  select(individual_id, timestamp, receiver_id, 
         receiver_easting, receiver_northing, receiver_range) |>
  as.data.table()

#### Align acoustic & archival time series
obs <- 
  lapply(split(acoustics_w_arc, acoustics_w_arc$individual_id), function(acc) {
    
    # acc <- split(acoustics_w_arc, acoustics_w_arc$individual_id)[["17"]]
    print(acc$individual_id[1])
    
    # Define archival data 
    arc <- 
      archival |> 
      filter(individual_id == acc$individual_id[1]) |> 
      filter(timestamp %within% lubridate::interval(min(acc$timestamp), max(acc$timestamp))) |> 
      arrange(timestamp) |> 
      as.data.table()
    
    # Fill gaps in archival time series
    # TO DO
    # * Evaluate behaviour of pf_forward() when depth = NA
    # * .update_ac() function needs to be able to handle this
    if (!all(unique(diff(arc$timestamp)) %in% 2)) {
      print("... Regularising archival time series...")
      fill <- data.table(individual_id = arc$individual_id[1], 
                         timestamp = seq(min(arc$timestamp), max(arc$timestamp), "2 mins"))
      arc <- left_join(fill, 
                       arc |> select(-individual_id), 
                       by = "timestamp")
    }
    
    # Collate observations
    obs <- acs_setup_obs(.acoustics = acc, .archival = arc, 
                         .step = "2 mins", 
                         .mobility = mobility, 
                         .detection_range = gamma)
    obs[, individual_id := acc$individual_id[1]]
    obs
    
  }) |> 
  rbindlist()

#### Identify detections & depths for analysis
# Identify detections 
detections <- 
  obs |> 
  filter(detection == 1L) |> 
  filter(!is.na(depth)) |> 
  select(individual_id, timestamp, detection_id, receiver_id, depth) |> 
  arrange(individual_id, timestamp) |>
  tidyr::unnest(cols = receiver_id) |>
  group_by(individual_id, detection_id) |> 
  arrange(receiver_id, .by_group = TRUE) |>
  mutate(receiver_key = paste(receiver_id, collapse = ", ")) |>
  mutate(receiver_id = list(receiver_id)) |>
  slice(1L) |>
  ungroup() |>
  as.data.table()
unique(detections$receiver_key)
nrow(detections)
# Assign depth limits
de <- calc_depth_error(detections$depth)
detections[, depth_shallow := detections$depth + de[1, ]]
detections[, depth_deep := detections$depth + de[2, ]]
head(detections)
# Check for duplicated detections 
detections |> 
  group_by(receiver_key) |> 
  summarise(dup = length(which(duplicated(depth)))) |> 
  arrange(dup) |> 
  as.data.table()
# Drop duplicated depth observations for speed
detsbt <- 
  detections |> 
  group_by(receiver_id) |> 
  filter(!duplicated(depth)) |> 
  as.data.table()
unique(detections$receiver_key)
length(unique(detections$receiver_key))

#### Count the number of possible locations in each detection container (~ 3 mins)
detsbt <- 
  pbapply::pblapply(split(detsbt, detsbt$receiver_key), cl = NULL, function(d) {
    
    # d <- split(detsbt, detsbt$receiver_key)[["26, 31"]]
    
    #### Define acoustic container 
    # Define receivers
    r <- 
      moorings |> 
      filter(receiver_id %in% d$receiver_id[[1]]) |> 
      as.data.table()
    # Define detection container (including the intersection)
    v <- 
      split(r, seq_row(r)) |>
      lapply(function(.r) {
        # .r <- r[1, ]
        cbind(.r$receiver_easting, .r$receiver_northing) |>
          terra::vect() |> 
          terra::buffer(width = .r$receiver_range + 250, quadsegs = 1e3)
      }) |> 
      patter:::spatIntersect()
    terra::crs(v) <- terra::crs(bathy)
    
    #### Expect depth values with container
    depth <- exactextractr::exact_extract(bathy, sf::st_as_sf(v))[[1]]$value
    
    #### Count the number of possible locations for each observation
    d$count <- 
      lapply(seq_row(d), function(i) {
        count <- table(depth >= d$depth_shallow[i] & depth <= d$depth_deep[i])["TRUE"]
        count <- as.numeric(count)
        if (is.na(count)) count <- 0
        count
      }) |> unlist()
    
    #### Return update data.table with counts
    d
    
  }) |> rbindlist()

#### Update detections
detections <- 
  detections |>
  left_join(detsbt |> select(receiver_key, depth, count), 
            by = c("receiver_key", "depth"))

#### Identify problematic receiver/depth combinations
# Problematic receivers: 30, (30 + 33), 39, 4, 40
detections[count <= 20L, ]
issues <- 
  detections |>
  select(individual_id, timestamp, receiver_key, depth, depth_shallow, depth_deep, count) |>
  filter(count < 20L) |> 
  arrange(receiver_key, depth, timestamp)

#### Check receiver locations
# Plot bathy 
terra::plot(bathy, 
            xlim = range(moorings$receiver_easting), 
            ylim = range(moorings$receiver_northing))
# Visualise the range of the Howe data 
terra::plot(howe, col = scales::alpha("grey", 0.5), add = TRUE)
# Add detection containers
cbind(moorings$receiver_easting, moorings$receiver_northing) |>
  terra::vect() |> 
  terra::buffer(width = 750) |> 
  terra::lines()
# Label receivers
basicPlotteR::addTextLabels(moorings$receiver_easting, 
                            moorings$receiver_northing, 
                            moorings$receiver_id)
# > receivers 30, 33 are in the southern receiver curtain
# > receiver 39 is south of Morvern (in Howe data)
# > receivers 4 & 40 are west lismore (in Digi data)

#### Compare observed depths to distribution of expected depths for problematic receiver(s)
# Define receiver coordinates
rid <- 40
r <- 
  moorings |> 
  filter(receiver_id %in% rid) |> 
  select(receiver_easting, receiver_northing) |> 
  as.matrix()
# Define detection container
v <- 
  lapply(seq_row(r), function(i) {
    r[i, , drop = FALSE] |> 
      terra::vect() |> 
      terra::buffer(width = gamma) 
  }) |> 
  patter:::spatIntersect()
# Record depths 
depths <- terra::extract(bathy, v)$depth
# Compare
range(issues$depth[issues$receiver_key == paste(rid, collapse = ", ")])
range(depths)

#### Results 
#
# * Gamma: 750 m:
#
# * Summary 
# > Individual(s) are 7 - 17 m shallower than the container (for receivers in Howe)
# > Individual(s) are 30 - 45 m shallower than the container (for receivers in Digi)
#
# * Receiver 30:
# > range(issues$depth[issues$receiver_key == paste(rid, collapse = ", ")]) # [OBS]
# [1] 11.56 16.40
# > range(depths) # [IN CONTAINER]
# [1]  23.71876 157.19899
# --> >= individual is 7 m shallower than container
#
# * Receiver 30 + 33:
# > range(issues$depth[issues$receiver_key == paste(rid, collapse = ", ")])
# [1] 12.99 12.99
# > range(depths)
# [1]  23.71876 150.09175
# --> >= 11 m shallower than container
#
# * Receiver 39: 
# > range(issues$depth[issues$receiver_key == paste(rid, collapse = ", ")])
# [1] 29.98 35.74
# > range(depths)
# [1]  52.72851 160.66391
# --> >= 17 m shallower than container
#
# * Receiver 4:
# > range(issues$depth[issues$receiver_key == paste(rid, collapse = ", ")])
# [1] 39.85 52.54
# > range(depths)
# [1]  97.57571 181.93729
# --> >= 45 m shallower than container
# 
# * Receiver 40:
# > range(issues$depth[issues$receiver_key == paste(rid, collapse = ", ")])
# [1] 69.68 69.68
# > range(depths)
# [1] 102.6166 185.1672
# --> >= 30 m shallower than container


###########################
###########################
#### Container dynamics

# (optional) TO DO
# * Validate there are possible locations within each acoustic container
# * For individuals that left the MPA, we expect some departure from this expectation


#### End of code. 
###########################
###########################