###########################
###########################
#### synthesis.R

#### Aims
# 1) Synthesis real-world analyses

#### Prerequisites
# 1) Previous scripts


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
# try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
dv::src()

#### Load data 
skateids          <- qs::qread(here_data("input", "mefs", "skateids.qs"))
unitsets          <- qs::qread(here_data("input", "unitsets.qs"))
moorings_in_mpa   <- qs::qread(here_data("input", "mefs", "moorings-in-mpa.qs"))
acoustics_raw     <- qs::qread(here_data("input", "mefs", "acoustics.qs"))
archival_raw      <- qs::qread(here_data("input", "mefs", "archival.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))
archival_by_unit  <- qs::qread(here_data("input", "archival_by_unit.qs"))
behaviour_by_unit <- qs::qread(here_data("input", "behaviour_by_unit.qs"))
recaps            <- readRDS(here_data_raw("movement", "recaptures_processed.rds"))
iteration_patter  <- qs::qread(here_data("input", "iteration", "patter.qs"))
mpa               <- qreadvect(here_data("spatial", "mpa.qs"))


###########################
###########################
#### Visualise time series 

# Check individual_id, acoustic_id, dst_id for comparison to Lavender et al. (2021)
# (Add sex, maturity etc. manually from Lavender et al. (2021))
skateids |>
  filter(individual_id %in% rbindlist(acoustics_by_unit)$individual_id) |> 
  mutate(acoustic_id = substr(acoustic_id, 4, 6)) |>
  select(individual_id, acoustic_id, dst_id) |>
  as.data.table() |> 
  prettyGraphics::tidy_write("./fig/individuals.txt")

if (FALSE) {
  
  #### Process modelled time series 
  # Define acoustics
  acoustics <- rbindlist(acoustics_by_unit)
  # Define archival dataset with behavioural state
  for (i in 1:length(archival_by_unit)) {
    if (!is.null(archival_by_unit[[i]])) {
      archival_by_unit[[i]][, behaviour := behaviour_by_unit[[i]]]
    }
  }
  archival <- rbindlist(archival_by_unit)
  archival[, depth := abs(depth) * -1]
  range(archival$depth)
  # Define keys for matching (below)
  acoustics[, key := paste(individual_id, mmyy)]
  archival[, key := paste(individual_id, mmyy)]
  # Collect time series
  timeseries <- bind_rows(archival, acoustics)
  
  #### Process 'raw' time series
  # Define 'raw' acoustics
  # * Focus on modelled individuals
  # * Focus on months within the study period 
  # * Focus on months _NOT_ included in models
  # * In the figure below, we add the raw datasets in grey 
  # * This helps to demonstrate the wider data context
  acoustics_raw <- 
    acoustics_raw |> 
    filter(individual_id %in% unitsets$individual_id) |> 
    mutate(mmyy = mmyy(timestamp), 
           key = paste(individual_id, mmyy)) |>
    filter(!(key %in% acoustics$key)) |> 
    filter(mmyy %in% unique(timeseries$mmyy)) |> 
    as.data.table()
  # Define 'raw' archival
  archival_raw <- 
    archival_raw |> 
    filter(individual_id %in% unitsets$individual_id) |> 
    mutate(depth = abs(depth) * -1, 
           mmyy = mmyy(timestamp), 
           key = paste(individual_id, mmyy)) |> 
    filter(!(key %in% archival$key)) |> 
    filter(mmyy %in% unique(timeseries$mmyy)) |> 
    as.data.table()
  # Define recaptures
  # * We mark recaptures as vertical lines on the plot
  # * This shows why some individuals with otherwise comphrensive time series have been excluded
  recaps <- 
    recaps |> 
    filter(code == "rc") |> 
    mutate(timestamp = as.POSIXct(paste(date, "12:00:00"), tz = "UTC")) |> 
    select(dst_id, timestamp) |> 
    mutate(individual_id = skateids$individual_id[match(dst_id, skateids$dst_id)]) |>
    filter(individual_id %in% timeseries$individual_id) |> 
    as.data.table()
  
  #### Define plot properties
  # x axis tick marks (days 5, 15 and 25 on each month)
  xticks <- seq(floor_date(min(timeseries$timestamp), "month"), 
                ceiling_date(max(timeseries$timestamp), "month"),
                by = "1 day")
  xticks <- xticks[day(xticks) %in% c(5, 15, 25)]
  
  #### Plot time series (~15 s)
  tic()
  png(here_fig("mefs", "time-series.png"), 
      height = 8, width = 10, res = 600, units = "in")
  ggplot(timeseries) +
    geom_line(data = archival_raw, aes(timestamp, depth), 
              lwd = 0.25, colour = "lightblue", alpha = 0.3) + 
    geom_text(data = acoustics_raw, 
              aes(timestamp, 2.5, label = "|"),
              colour = "grey70", size = 2, vjust = 0.5) + 
    geom_line(data = archival, 
              aes(timestamp, depth, colour = depth), 
              lwd = 0.25) + 
    geom_text(data = acoustics,
              aes(timestamp, 2.5, label = "|"),
              colour = "black", size = 2, vjust = 0.5) + 
    geom_vline(data = recaps, aes(xintercept = timestamp), 
               colour = "red", alpha = 0.5,  linewidth = 0.35) + 
    scale_x_datetime(labels = scales::date_format("%d"), 
                     breaks = xticks, 
                     expand = c(0.025, 0)) +
    scale_y_continuous(breaks = c(-50, -250), expand = c(0, 0), limits = c(-300, 5)) + 
    scale_color_gradientn(
      colours = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(256)),
      limits = c(-350, 0),
      na.value = "grey"
    ) + 
    xlab("Time (day of month)") + 
    ylab("Depth (m)") + 
    facet_grid(individual_id ~ mmyy, scales = "free_x") + 
    theme_bw() + 
    theme(legend.position = "none", 
          panel.grid = element_blank(),
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(r = 10)))
  dev.off()
  toc()
 
  #### Selected checks
  acoustics |> 
    filter(individual_id == 35 & mmyy == "11-2016") |> 
    print(125)
  
  #### Summary statistics
  # Number of detections per month
  ndet <- 
    acoustics |> 
    group_by(individual_id, mmyy) |> 
    summarise(n = n()) |> 
    arrange(n) |>
    ungroup() |>
    as.data.table()
  median(ndet$n)
  # Gaps between sequential detections
  gdet <- 
    acoustics |> 
    group_by(individual_id, mmyy) |> 
    mutate(gap = serial_difference(timestamp, units = "mins")) |> 
    summarise(utils.add::basic_stats(gap, na.rm = TRUE)) |> 
    ungroup() |>
    as.data.table()
  gdet
  max(gdet$max)
}


###########################
###########################
#### Patterns of space use

#### (Quick) Plot MPA
terra::plot(mpa, xlim = terra::ext(mpa)[1:2], ylim = terra::ext(mpa)[3:4])
terra::sbar(10000)

#### Example UDs from each algorithm
algorithms <- c("COA", "RSP", "AC", "DC", "ACDC")
mapfiles <-
  data.table(algorithm = factor(algorithms, algorithms),
             mapindex  = seq_len(length(algorithms)),
             mapfile   = c("coa/ac/1/ud/spatstat/h/ud.tif", 
                           "rsp/ac/1/ud/dbbmm/ud.tif", 
                           "patter/ac/1/ud/spatstat/h/ud.tif", 
                           "patter/dc/1/ud/spatstat/h/ud.tif", 
                           "patter/acdc/1/ud/spatstat/h/ud.tif"))
mapfiles <- 
  CJ(individual_id = c(25, 35, 36), # c(13, 20, 25, 27, 29, 35, 36, 38),
     mapindex      = mapfiles$mapindex
  ) |>  
  mutate(row    = factor(individual_id, levels = sort(unique(individual_id))), 
         column = mapfiles$algorithm[match(mapindex, mapfiles$mapindex)],
         mapfile = mapfiles$mapfile[match(mapindex, mapfiles$mapindex)],
         mapfile = file.path("data", "output", "analysis", individual_id, "04-2016", mapfile)) |> 
  select(row, column, mapfile) |> 
  as.data.table()
ggplot_maps(mapfiles, 
            png_args = list(filename = here_fig("analysis", "ud-examples.png"), 
                            height = 4, width = 4, units = "in", res = 800))

#### Overall ACDC UD with tagging locations & angling records
# List tif files
mapfiles <- 
  iteration_patter |> 
  filter(dataset == "acdc" & sensitivity == "best") |> 
  mutate(mapfile = file.path(folder_ud, "spatstat", "h", "ud.tif")) |> 
  as.data.table()
# Read UDs
uds <- 
  lapply(mapfiles$mapfile, function(ud.tif) {
  if (file.exists(ud.tif)) {
    terra::rast(ud.tif)
  }
}) |> 
  plyr::compact()
# Compute 'overall' UD
uds <- do.call(c, uds)
ud <- terra::app(uds, fun = "sum", na.rm = TRUE) / terra::nlyr(uds)
stopifnot(all.equal(1, as.numeric(terra::global(ud, "sum", na.rm = TRUE))))
ud.tif <- tempfile(fileext = ".tif")
terra::writeRaster(ud, ud.tif)
# Get capture locations for relevant individuals 
tagsf <- 
  skateids |> 
  filter(individual_id %in% mapfiles$individual_id)|> 
  select(lon = long_tag_capture, lat = lat_tag_capture) |>
  as.matrix() |> 
  terra::project(from = "EPSG:4326", to = terra::crs(mpa)) |>
  as.data.frame() |> 
  sf::st_as_sf(coords = c("V1", "V2"), crs = terra::crs(mpa))
  # Get all angling records (download from 28/11/2024)
# (We expect ! NAs introduced by coercion warning here)
cr <- 
  here_data_raw("skatespotter", "data (15).csv") |> 
  read.csv() |> 
  select(individual_id, date = date_captured, lon = longitude, lat = latitude) |> 
  mutate(date = as.Date(date), 
         lon = as.numeric(lon), 
         lat = as.numeric(lat)) |>
  filter(!is.na(lon) & !is.na(lat)) |> 
  filter(lon != 0 & lat != 0) |>
  as.data.table()
range(cr$date)
p <- 
  ggplot(data.frame(year = lubridate::year(cr$date))) + 
  geom_histogram(aes(year), bins = 50) + 
  theme_bw() 
plotly::ggplotly(p)
# Get angling locations
crxy <- 
  cbind(cr$lon, cr$lat) |> 
  terra::vect(crs = "EPSG:4326") |> 
  terra::project(terra::crs(mpa)) |> 
  terra::crds(df = TRUE)
crsf <- 
  sf::st_as_sf(crxy, coords = c("x", "y"), crs = terra::crs(mpa))
# Quick plot
terra::plot(ud)
points(crxy$x, crxy$y)
hr <- map_hr_home(ud, .add = TRUE)
poly <- terra::as.polygons(hr == 1)
poly <- poly[poly[[1]] == 1]
poly <- poly |> sf::st_as_sf()
# ggplot
png(here_fig("analysis", "ud-overall.png"), 
    height = 3.5, width = 2.5, units = "in", res = 800)
p <- 
  ggplot_maps(data.table(mapfile = ud.tif, row = 1, column = 1), 
              png_args = NULL) 
  p + 
  geom_sf(data = poly, fill = NA) + 
  geom_sf(data = crsf, shape = 21, size = 0.001, linewidth = 0, colour = "purple", alpha = 0.2) + 
  geom_sf(data = tagsf, shape = 11, size = 0.9, colour = "darkgreen") + 
  coord_sf(xlim = p$coordinates$limits$x, ylim = p$coordinates$limits$y) + 
  theme(panel.grid.major = element_line(colour = "#0000FF", linewidth = 0.025))
dev.off()

#### Algorithm sensitivity


###########################
###########################
#### Residency 

#### Null model
res_null <- 
  qs::qread(here_data("output", "simulation-summary", "residency-null.qs")) |>
  slice(1:3) |> 
  select("zone", "time") |>
  mutate(zone = factor(zone, levels = c("open", "closed", "total"), labels = c("Open", "Closed", "Protected"))) |>
  as.data.table()
res_null

#### Detection days
# (Code modified from simulate-algorithms.R)
residency <- lapply(acoustics_by_unit, function(acoustics) {
  
  if (is.null(acoustics)) {
    return(NULL)
  }
  
  # Total number of days in month 
  # acoustics <- acoustics_by_unit[[1]]
  ndays <- as.integer(lubridate::days_in_month(as.Date.mmyy(acoustics$mmyy[1])))
  
  # Compute detection days for receivers in MPA
  dds_total <- 
    acoustics |>
    filter(receiver_id %in% moorings_in_mpa$receiver_id) |> 
    mutate(day = lubridate::day(timestamp)) |> 
    summarise(time = length(unique(day)) / ndays) |>
    mutate(individual_id = acoustics$individual_id[1], 
           month_id      = acoustics$mmyy[1], 
           unit_id       = acoustics$unit_id[1], 
           algorithm     = "DD", 
           sensitivity   = "Best",
           zone          = "total") |> 
    select(individual_id, month_id, unit_id, algorithm, sensitivity, zone, time) |> 
    arrange(individual_id, month_id, unit_id, algorithm, sensitivity, zone) |>
    as.data.table()
  
  # Detection days in closed areas are identical 
  # * (All receivers were in closed areas)
  # * For the figures, we also record DDs in closed areas
  dds_closed <- copy(dds_total)
  dds_closed[, zone := "closed"]
  
  rbind(dds_total, dds_closed)
  
}) |> rbindlist()

qs::qsave(residency, 
          here_data("output", "analysis-summary", "residency-detection-days.qs"))

#### Collate residency estimates
residency <- 
  rbindlist(
    list(
      qs::qread(here_data("output", "analysis-summary", "residency-detection-days.qs")),
      qs::qread(here_data("output", "analysis-summary", "residency-coa.qs")),
      qs::qread(here_data("output", "analysis-summary", "residency-rsp.qs")),
      qs::qread(here_data("output", "analysis-summary", "residency-patter.qs"))
    )
  ) |> 
  filter(!is.na(zone)) |>
  mutate(month = as.Date.mmyy(month_id)) |>
  mutate(individual_id = factor(individual_id)) |>
  mutate(zone = factor(zone, levels = c("open", "closed", "total"), labels = c("Open", "Closed", "Protected"))) |>
  as.data.table()

#### Compute summary statistics
residency |> 
  group_by(zone, algorithm) |>
  summarise(min = min(time), 
            max = max(time)) |>
  ungroup() |>
  arrange(zone, algorithm, min)
residency |> 
  filter(algorithm == "DD") |> 
  summarise(utils.add::basic_stats(time, na.rm = TRUE))
residency |> 
  filter(algorithm %in% c("COA", "RSP")) |> 
  filter(sensitivity == "Best") |>
  group_by(zone) |>
  summarise(utils.add::basic_stats(time, na.rm = TRUE))
residency |> 
  filter(algorithm == "AC") |> 
  filter(sensitivity == "Best") |>
  group_by(zone) |>
  summarise(utils.add::basic_stats(time, na.rm = TRUE))
residency |> 
  filter(algorithm == "DC") |> 
  filter(sensitivity == "Best") |>
  group_by(zone) |>
  summarise(utils.add::basic_stats(time, na.rm = TRUE))
residency |> 
  filter(algorithm == "ACDC") |> 
  # filter(sensitivity == "Best") |>
  group_by(zone) |>
  summarise(utils.add::basic_stats(time, na.rm = TRUE))

#### Visualise residency trends
head(residency)
png(here_fig("analysis", "residency-best.png"), 
    height = 6, width = 10, units = "in", res = 600)
residency |>
  # filter(sensitivity == "Best") |>
  # filter(zone == "total") |> 
  as_tibble() |> 
  ggplot() + 
  geom_point(aes(
    month, time, 
    colour = individual_id, 
    shape = sensitivity, 
    alpha = if_else(sensitivity == "Best", 1, 0.5), 
    size = if_else(sensitivity == "Best", 1, 0.5)
  )) + 
  geom_line(aes(
    month, time, 
    colour = individual_id, 
    group = interaction(individual_id, sensitivity), 
    alpha = if_else(sensitivity == "Best", 1, 0.5), 
    size = if_else(sensitivity == "Best", 0.75, 0.25)
  )) +
  scale_alpha_identity() +
  scale_size_identity() +
  scale_shape_manual(values = c(20, 17, 15, 3, 4, 8, 13)) + 
  scale_x_date(labels = scales::date_format("%b-%y")) +
  scale_y_continuous(expand = c(0, 0.05), limits = c(0, 1)) + 
  geom_hline(data = res_null, aes(yintercept = time, colour = NULL), linetype = "dashed") +
  xlab("Time (months)") + 
  ylab("Residency") + 
  guides(
    colour = guide_legend(order = 1, title = "Individual"),
    shape = guide_legend(order = 2, title = "Parameterisation"),
    alpha = "none",
    size = "none"
  ) + 
  facet_grid(zone ~ algorithm) + 
  theme_bw() + 
  theme(panel.spacing.y = unit(1.5, "lines"), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(), 
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)), 
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


#### End of code. 
###########################
###########################
