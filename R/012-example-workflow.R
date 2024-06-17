###########################
###########################
#### example-workflow.R

#### Aims
# 1) Prepare flapper algorithm inputs

#### Prerequisites
# 1) Obtain raw data
# 2) https://github.com/edwardlavender/flapper_appl/blob/master/R/define_global_param.R


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()
# op <- options(error = recover)
op <- options(error = function(...) beepr::beep(7))

#### Essential packages
library(dv)
library(data.table)
library(dtplyr)
library(dplyr)
library(arrow)
library(ggplot2)
library(patter)
library(prettyGraphics)
library(sf)
library(testthat)
library(tictoc)

#### Load data
dv::src()
bathy     <- terra::rast(here_data("spatial", "bathy.tif")) 
bset      <- terra::rast(here_data("spatial", "bset.tif")) 
im        <- qs::qread(here_data("spatial", "im.qs"))
win       <- qs::qread(here_data("spatial", "win.qs"))
moorings  <- readRDS(here_data("mefs", "moorings.rds"))
acoustics <- readRDS(here_data("mefs", "acoustics.rds"))
archival  <- readRDS(here_data("mefs", "archival.rds"))
pars      <- readRDS(here_data("input", "pars.rds"))
overlaps  <- readRDS(here_data("input", "overlaps.rds"))
kernels   <- acs_setup_detection_kernels_read()
ewin      <- readRasterLs(here_data("input", "depth-window"), index = FALSE)
coast     <- readRDS(here_data("spatial", "coast.rds")) |> st_geometry()

#### Local pars
seed <- 1L
run        <- FALSE
run_origin <- FALSE
rerun      <- FALSE
manual     <- run


###########################
###########################
#### Set up analysis

#### Process SpatRasters
terra:::readAll(bathy)
terra:::readAll(bset)

#### Process movement datasets
# Limits
xlim <- as.POSIXct(c("2016-07-01 00:00:00", "2016-08-01 23:58:00"), 
                   tz = "UTC")
# Acoustics
acc <- 
  acoustics |>
  filter(individual_id == 25) |> 
  filter(timestamp >= xlim[1]) |> 
  filter(timestamp <= xlim[2]) |> 
  as.data.table()
# Archival data (if applicable)
arc <- 
  archival |>
  filter(individual_id == 25) |> 
  filter(timestamp >= xlim[1]) |> 
  filter(timestamp <= xlim[2]) |> 
  arrange(timestamp) |>
  # Define 'resting': 0L = resting; 1L = not resting;
  mutate(va = abs(serial_difference(depth)), 
         state = if_else(va <= 0.25, true = 0L, false = 1L)) |>
  as.data.table()
arc$va <- NULL
arc$state[nrow(arc)] <- 1L

#### Build data list
dlist <- pat_setup_data(.acoustics = acc, 
                        .moorings = moorings, 
                        .archival = arc, 
                        .bathy = bathy, 
                        .lonlat = FALSE)
# Include AC likelihood terms (required for ACPF, ACDCPF)
dlist$algorithm$detection_overlaps <- overlaps
dlist$algorithm$detection_kernels  <- kernels
#
# TO DO
# * Time trial approaches bring kernels into memory
# * Identify the receivers with detections & overlapping receivers
# * Optionally crop & use xy coordinates to querying
#
# Include DC likelihood terms (required for DCPF, ACDCPF)
dlist$spatial$bset      <- bset
dlist$algorithm$ewindow <- ewin
dlist$algorithm$n       <- 1e3L
# Include parameters
dlist$pars$shape     <- pars$patter$shape
dlist$pars$scale     <- pars$patter$scale
dlist$pars$mobility  <- pars$patter$mobility
dlist$pars$gamma     <- pars$patter$detection_range
dlist$algorithm$dlen <- dtruncgamma
# Include additional elements below
# * .$spatial$origin element 
# * .$algorithm$pos_detections element (for acs_filter_container_acdc())

#### Define observations timeline
obs <- pf_setup_obs(.dlist = dlist,
                    .trim = FALSE,
                    .step = "2 mins", 
                    .mobility = 500, 
                    .receiver_range = 750)
# Include behavioural states
obs[, state := arc$state[match(timestamp, arc$timestamp)]]
# Include elements for acs_setup_container_acdc()
# * Define receiver_id_next_key (container) as in acs_setup_containers_rcd()
obs[, receiver_id_next_key := 
      lapply(obs$receiver_id_next, acs_setup_receiver_key) |> unlist()]
# * Define pos_detections
dlist$algorithm$pos_detections <- which(!sapply(obs$receiver_id, is.null))

#### Examine observations
# Check obs
obs
range(obs$timestamp)
stopifnot(any(obs$detection == 1L))
stopifnot(!any(is.na(obs$depth)))
# Visualise obs
plot(arc$timestamp, arc$depth*-1, 
     xlim = xlim, ylim = c(-225, 0),
     type = "n", xlab = "Time (month/day)", ylab = "Depth (m)")
s <- seq(nrow(arc) - 1)
arrows(arc$timestamp[s], arc$depth[s] * -1, 
       arc$timestamp[s + 1L], arc$depth[s + 1] * -1, 
       col = factor(arc$state), length = 0, lwd = 0.5)
points(acc$timestamp, rep(0L, nrow(acc)), col = "red")

#### Define origin (~40 s)
# This is included in dlist below, if necessary
origin <- NULL
if (TRUE) {
  tic()
  obs$depth[1]
  stopifnot(!is.na(obs$depth[1]))
  origin <- dc_origin(.ewindow = dlist$algorithm$ewindow, .depth = obs$depth[1])
  # terra::plot(origin)
  # terra::global(terra::mask(bathy, origin), "range", na.rm = TRUE)
  toc()
}

#### Define algorithm 
algs <- c("acpf", "acdcpf", "dcpf")
alg <- algs[2]


###########################
###########################
#### Forward filter


###########################
#### Set up 

#### Define directories
pff_folder <- here_data("example", "forward", alg, "output")
log.txt    <- here_data("example", "forward", alg, "log.txt")
# Optionally wipe old directories
if (run) {
  unlink(log.txt)
  unlink(pff_folder, recursive = TRUE)
}
# (Re)build
if (!dir.exists(pff_folder)) {
  dir.create(pff_folder, recursive = TRUE)
}

#### Define data list
if (alg %in% c("dcpf", "acdcpf")) {
  dlist$spatial$origin <- origin
} else {
  dlist$spatial$origin <- NULL
}

#### Define args
# Define movement args 
# * TO DO
# * Define behaviourally dependent movement model (via .rpropose and .dpropose)
# * Downgrade movement jumps using shortest distances model
if (pars$patter$model == "truncated gamma") {
  margs <- list(.shape = dlist$pars$shape, 
                .scale = dlist$pars$scale, 
                .mobility = dlist$pars$mobility)
  rargs <- margs
  dargs <- margs
} else if (pars$patter$model == "uniform") {
  rargs <- list(.rlen = rlenuniform, .mobility = dlist$pars$mobility)
  dargs <- list(.dlen = dlenuniform, .mobility = dlist$pars$mobility)
} else {
  stop("Input to `pars$patter$model` not supported.")
}
# Define trial arguments
trial <- 
  pf_opt_trial(
    .trial_kick = 2L, 
    .trial_sampler_crit = 50L, 
    .trial_revert = 0L
  )
# Define record options
# c("timestep", "cell_past", "cell_now", "x_now", "y_now", "lik")
record <- 
  pf_opt_record(
    .save = FALSE,
    .sink = pff_folder, 
    .cols = NULL 
  )
# Collate arguments
args <- list(.obs = obs,
             .dlist = dlist,
             .rpropose = pf_rpropose_kick, 
             .dpropose = pf_dpropose,
             .rargs = rargs, .dargs = dargs,
             .n = dlist$algorithm$n, 
             .trial = trial,
             .record = record,
             .control = pf_opt_control(.sampler_batch_size = 1000L),
             .verbose = log.txt) 
# Define algorithm-specific likelihoods
if (alg == "acpf") {
  args$.likelihood <- list(pf_lik_ac = pf_lik_ac, 
                           acs_filter_land = acs_filter_land, 
                           acs_filter_container = acs_filter_container)
} else if (alg == "dcpf") {
  args$.likelihood <- list(pf_lik_dc = pf_lik_dc_2)
} else if (alg == "acdcpf") {
  args$.likelihood <- list(pf_lik_dc = pf_lik_dc_2, 
                           pf_lik_ac = pf_lik_ac, 
                           acs_filter_container = acs_filter_container)
                           # acs_filter_container = acs_filter_container_acdc)
}


###########################
#### Forward run 

#### Timings
# * ACPF: ~107 mins
# * DCPF: 
# * ACDCPF: 

#### Initial run 
if (run) {
  tic()
  set.seed(seed)
  out_pff <- do.call(patter::pf_forward, args)
  toc()
  beepr::beep(10L)
  saveRDS(out_pff, file.path(pff_folder, "out_pff.rds"))
} else {
  out_pff <- readRDS(file.path(pff_folder, "out_pff.rds"))
}

#### (optional) Rerun
if (rerun) {
  rerun_args <- args
  rerun_args$.rerun      <- out_pff
  rerun_args$.rerun_from <- 12000
  # nt <- length(pf_files(file.path(pff_folder, "history")))
  # rerun_args$.rerun_from <- plyr::round_any(nt - 500, 500, floor)
  tic()
  set.seed(seed)
  out_pff_2 <- do.call(patter::pf_forward, rerun_args)
  toc()
  beepr::beep(10L)
}

# To debug convergence issues, see ./R/supporting/convergence/.
# * 4572     
# * 9267/9270
# * 12510     - can be resolved by forward/backward reversal

# Warning: Convergence error: there are no particles with positive weights at time step 12493. Returning outputs up to time step 12493.


###########################
#### Outputs

#### Output file size (MB)
hs <- pf_files_size(pff_folder, .folder = "history")
ds <- pf_files_size(pff_folder, .folder = "diagnostics")
sum(c(hs, ds))

#### Output diagnostics 
if (manual) {
  diag <- pf_diag_convergence(pff_folder)
  tail(diag, 100)
}


###########################
#### Validation 

if (manual) {
  
  #### Collate particle samples
  h <- patter:::.pf_history_dt(file.path(pff_folder, "history"), 
                               schema = arrow::schema(timestep = int32(), 
                                                      cell_past = int32(), 
                                                      x_past = double(), 
                                                      y_past = double(), 
                                                      cell_now = int32(), 
                                                      x_now = double(),
                                                      y_now = double(), 
                                                      lik = double()))
  head(h)
  
  #### Validate locations
  # Validate grid cells & coordinates align
  expect_equal(
    as.matrix(h[, .(x_now, y_now)]) |> unname(),
    terra::xyFromCell(bathy, h$cell_now) |> unname())
  
  #### Validate proposals (distances)
  # Validate distances are < mobility (accounting for grid resolution)
  xy0 <- as.matrix(h[, .(x_now, y_now)])
  xy1 <- terra::xyFromCell(bathy, h$cell_past)
  len <- terra::distance(xy0, xy1, lonlat = FALSE, pairwise = TRUE)
  range(len, na.rm = TRUE)
  expect_true(all(len < 500 + terra::res(bathy)[1]/2, na.rm = TRUE))
  
  #### Validate detection information 
  # At the moment of detection, particle samples should be within detection ranges
  hd <- h[timestep %in% obs$timestep[obs$detection == 1L]]
  hd <- 
    hd |> 
    left_join(
      obs |> 
        select(timestep, receiver_id),
      by = "timestep") |> 
    tidyr::unnest(cols = receiver_id) |> 
    left_join(dlist$data$moorings |> 
                select(receiver_id, receiver_x, receiver_y), 
              by = c("receiver_id")) |> 
    mutate(dist = terra::distance(cbind(x_now, y_now), 
                                  cbind(receiver_x, receiver_y), 
                                  lonlat = FALSE, 
                                  pairwise = TRUE)) |> 
    as.data.table()
  range(hd$dist)
  expect_true(all(hd$dist < 750))
  # In detection gaps, particle samples should (generally) be beyond detection ranges
  # * TO DO
  
  #### Validate depth information 
  h$bathy <- terra::extract(bathy, h$cell_now)
  h$depth <- obs$depth[match(h$timestep, obs$timestep)]
  hist(h$bathy - h$depth)
  range(h$bathy - h$depth)
  
  #### Validate algorithm parameters
  # The number of particles at each time step
  np <- 
    h |>
    count(timestep) |>
    pull(n)
  expect_true(all(np == 1000L))
  
}


###########################
###########################
#### Backward sampler

#### Define directories
pfb_folder  <- here_data("example", "backward", alg)
pfbk_folder <- file.path(pfb_folder, "killer", "output")
log.txt     <- file.path(pfb_folder, "killer", "log.txt")
if (FALSE) {
  unlink(pfbk_folder, recursive = TRUE)
  unlink(log.txt)
}
if (!dir.exists(pfbk_folder)) {
  dir.create(pfbk_folder, recursive = TRUE)
}
record$sink <- pfbk_folder

#### Implement backward killer 
# ACPF: ~140 s
if (run) {
  tic()
  out_pfbk <- pf_backward_killer(.history = file.path(pff_folder, "history"),
                                 .record = record,
                                 .verbose = log.txt)
  toc()
}

#### Diagnostics 
if (manual) {
  # Extract diagnostics 
  # * 1 CPU: ~38 s
  # * 10 CPU (fork): 7 s (with chunks)
  pfbk_diag <- pf_diag_summary(pfbk_folder, 
                               schema = schema(timestep = int32(), 
                                               cell_now = int32(), 
                                               weight = double()))
  # Examine diagnostic traces
  plot(pfbk_diag$timestep, pfbk_diag$n_u, type = "l")
  plot(pfbk_diag$timestep, pfbk_diag$n_u, ylim  = c(0, 100), type = "l")
  plot(pfbk_diag$timestep, pfbk_diag$ess, type = "l")
}


###########################
###########################
#### Backward sampler

#### Define directories 
pfbs_folder <- file.path(pfb_folder, "sampler", "output")
log.txt     <- file.path(pfb_folder, "sampler", "log.txt")

#### Run sampler 
files <- pf_files(file.path(pff_folder, "history"))
files <- files[(length(files) - 10L):length(files)]
if (run) {
  tic()
  set.seed(seed)
  # TO DO: update!
  out_pfbs <- pf_backward_sampler(files, 
                                  .step_dens = dstep, lonlat = FALSE,
                                  .write_history = list(sink = pfbs_folder), 
                                  .txt = log.txt)
  toc()
}


###########################
###########################
#### Particle histories

#### Define plot region (~8 s)
# pxy <- pf_coord(file.path(pff_folder, "history"), .bathy = bathy)
tic()
ext <- 
  pff_folder |> 
  file.path("history") |> 
  arrow::open_dataset() |> 
  summarise(xmin = min(x_now), xmax = max(x_now), ymin = min(y_now), ymax = max(y_now)) |> 
  as.data.table() |> 
  summarise(xmin = min(xmin), xmax = max(xmax), ymin = min(ymin), ymax = max(ymax)) |>
  as.data.table()
ext <- terra::ext(ext$xmin, ext$xmax, ext$ymin, ext$ymax)
# Inflat ext
ext <- ext + 1e3
toc()

# Crop grid accordingly 
dlist_cpy <- dlist
dlist_cpy$spatial$bathy <- terra::crop(dlist_cpy$spatial$bathy, ext)
terra::plot(dlist_cpy$spatial$bathy)
points(moorings$receiver_easting, moorings$receiver_northing)

#### Define particle samples from the forward and backward runs
ff <- pf_files(pff_folder, "history")
fb <- pf_files(pfbk_folder)

#### Plot parameters
# Number of time steps
nt        <- length(ff)
# Observations 
obsp      <- obs[seq_len(nt), ]
obsp[, depth_neg := depth * -1]
# Receivers that recorded detections
receivers <- unique(unlist(obsp$receiver_id))
# Corresponding detection containers
v <- 
  moorings |> 
  filter(receiver_id %in% receivers) |> 
  dplyr::select(receiver_easting, receiver_northing) |> 
  as.matrix() |> 
  terra::vect() |> 
  terra::buffer(width = pars$patter$detection_range)
# Define receiver clumps
clumps <- 
  obsp |>
  dplyr::select(timestamp, detection, receiver_id) |> 
  filter(detection == 1L) |> 
  tidyr::unnest(cols = c(receiver_id)) |> 
  mutate(clump = rleid(receiver_id)) |>
  group_by(clump) |> 
  slice(1L) |> 
  ungroup() |> 
  as.data.table()
# Bathymetry colours
zlim      <- c(0, 250) 
zlim_neg  <- sort(zlim * -1)
col_param <- pretty_cols_brewer(zlim = zlim, 
                                scheme = "Blues")
col_df <- data.frame(from = col_param$breaks[1:(length(col_param$breaks) - 1)], 
                     to = col_param$breaks[2:(length(col_param$breaks))], 
                     color = col_param$col)

#### Make maps (~8 mins, 10 cl; ~1 hr 16, 0 cl)
if (FALSE) {
  png_sink <- here_fig("example", "animation", "frames", "maps")
  dir.create(png_sink, recursive = TRUE)
  tic()
  pf_plot_history(.dlist = dlist_cpy, 
                  .forward = ff, .backward = fb, 
                  .steps = seq_len(nt), 
                  .add_surface = list(range = col_param$zlim, col = col_df, legend = FALSE),
                  .add_forward = list(pch = "."), 
                  .add_layer = function() {
                    # Add detection curtains around relevant receivers
                    terra::lines(v)
                    # Add coast & moorings 
                    plot(coast, add = TRUE, col = scales::alpha("palegreen3", 0.5))
                    text(moorings$receiver_easting, moorings$receiver_northing, moorings$receiver_id, 
                         font = 2)
                  }, 
                  .png = list(filename = png_sink, height = 7, width = 7, units = "in", res = 300), 
                  .cl = 10L)
  toc()
}

#### Make time series (2 mins s, 10 cl)
if (FALSE) {
  png_sink <- here_fig("example", "animation", "frames", "obs")
  dir.create(png_sink, recursive = TRUE)
  tic()
  cl_lapply(seq_len(nt), .cl = 10L, .fun = function(t) {
    # Blank plot 
    png(file.path(png_sink, glue::glue("{t}.png")), 
        height = 7, width = 7, units = "in", res = 300)
    pretty_plot(obsp$timestamp, obsp$depth_neg, 
                xlab = "", ylab = "",
                pretty_axis_args = list(side = 3:2, 
                                        axis = list(list(format = "%d-%b"), 
                                                    list())),
                type = "n")
    # Add depth time series
    add_lines(x = obsp$timestamp, 
              y1 = obsp$depth_neg, 
              y2 = obsp$depth,
              breaks = col_param$breaks, cols = col_param$col, 
              lwd = 2)
    # Add detections
    p <- which(obsp$detection == 1L)
    points(obsp$timestamp[p], rep(0, length(p)))
    text(clumps$timestamp, rep(-20, nrow(clumps)), clumps$receiver_id, font = )
    # basicPlotteR::addTextLabels(as.numeric(obsp$timestamp[p]), rep(0, length(p)), labels = obsp$receiver_id[p] |> unlist())
    # Add time step
    lines(rep(obsp$timestamp[t], 2), zlim_neg, col = "red", lwd = 1.5)
    dev.off()
    NULL
  })
  toc()
}

#### Make videos(s)
if (FALSE) {
  
  #### Input files
  png_sink   <- here_fig("example", "animation", "frames", "maps")
  input      <- unlist(pf_files(png_sink))
  # input      <- input[c(1:500)]
  
  #### Parameters
  # Framerate
  fr <- 100
  # Total video duration:
  length(input) / fr
  
  #### GIF
  # Gifs only work well with delay > 0.005 s, which is too slow:
  # output     <- file.path(dirname(png_sink), "ani-1.gif")
  # gifski::gifski(input, output, delay = 0.001, width = 672, height = 672)
  
  #### Map mp4 (~9.3 mins)
  tic()
  output     <- here_fig("example", "animation", "ani-1.mp4")
  av::av_encode_video(input, output, framerate = fr)
  toc()

  #### Time series mp4 (~5 mins)
  tic()
  png_sink   <- here_fig("example", "animation", "frames", "obs")
  input      <- unlist(pf_files(png_sink))[seq_len(length(input))]
  output     <- here_fig("example", "animation", "ani-2.mp4")
  av::av_encode_video(input, output, framerate = fr)
  toc()

}


###########################
###########################
#### Maps

#### Define input folder
folder <- pfbk_folder
# folder <- pfbs_folder

#### Collect coordinates (~5 s)
# * We collect coordinates and use pf_map_dens()
# * This is faster than pf_map_pou() plus pf_map_dens()
tic()
pxy <- pf_coords(folder, .bathy = bathy)
toc()
if (manual) {
  # Number of unique cells sampled
  length(unique(pxy$cell_id))
  # Unique samples through time
  pxy |> 
    group_by(timestep) |> 
    summarise(n = n()) |> 
    ggplot() + 
    geom_line(aes(timestep, n))
}

#### Density (~70 s)
tic()
dens <- pf_map_dens(.xpf =  bathy, 
                    .im = im, 
                    .owin = win, 
                    .coord = pxy
                    )
if (manual) {
  # Add unique particle samples
  pxyu <- pxy[!duplicated(cell_id), ]
  points(pxyu$cell_x, pxyu$cell_y, 
         pch = ".", 
         col = scales::alpha("darkgrey", 0.7))
}
toc()

#### Write maps to file (~11 s)
tic()
terra::writeRaster(dens, 
                   here_data("example", "map", "acpf", "dens.tif"), 
                   overwrite = TRUE)
toc()


#### End of code
###########################
###########################