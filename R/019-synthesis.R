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
Sys.setenv("JULIA_SESSION" = "FALSE")
dv::src()

#### Load data 
bathy             <- terra::rast(here_data("spatial", "bathy.tif"))
bathy_5m          <- terra::rast(here_data("spatial", "bathy-5m.tif"))
mpa               <- qreadvect(here_data("spatial", "mpa.qs"))
skateids          <- qs::qread(here_data("input", "mefs", "skateids.qs"))
recaps            <- readRDS(here_data_raw("movement", "recaptures_processed.rds"))
moorings          <- qs::qread(here_data("input", "mefs", "moorings.qs"))
moorings_in_mpa   <- qs::qread(here_data("input", "mefs", "moorings-in-mpa.qs"))
acoustics_raw     <- qs::qread(here_data("input", "mefs", "acoustics.qs"))
archival_raw      <- qs::qread(here_data("input", "mefs", "archival.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))
archival_by_unit  <- qs::qread(here_data("input", "archival_by_unit.qs"))
behaviour_by_unit <- qs::qread(here_data("input", "behaviour_by_unit.qs"))
unitsets          <- qs::qread(here_data("input", "unitsets.qs"))
iteration_coa     <- qs::qread(here_data("input", "iteration", "coa.qs"))
iteration_rsp     <- qs::qread(here_data("input", "iteration", "rsp.qs"))
iteration_patter  <- qs::qread(here_data("input", "iteration", "patter-tidy.qs"))


###########################
###########################
#### Visualise time series 

if (FALSE) {
  
  # Check individual_id, acoustic_id, dst_id for comparison to Lavender et al. (2021)
  # (Add sex, maturity etc. manually from Lavender et al. (2021))
  skateids |>
    filter(individual_id %in% rbindlist(acoustics_by_unit)$individual_id) |> 
    mutate(acoustic_id = substr(acoustic_id, 4, 6)) |>
    select(individual_id, acoustic_id, dst_id) |>
    as.data.table() |> 
    prettyGraphics::tidy_write("./fig/individuals.txt")
  
  
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
  png(here_fig("mefs", "time-series-incl-legend.png"), 
      height = 8, width = 10, res = 1200, units = "in")
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
      breaks = seq(-350, 0, length.out = 100),
      colours = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(256)),
      limits = c(-350, 0),
      na.value = "grey"
    ) + 
    xlab("Time (day of month)") + 
    ylab("Depth (m)") + 
    facet_grid(individual_id ~ mmyy, scales = "free_x") + 
    theme_bw() + 
    theme(# legend.position = "none", 
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

###########################
#### Study area

#### (Quick) Plot study area
if (FALSE) {
  terra::plot(terra::rast(here_data("spatial", "bathy-5m.tif")))
  qs::qread(here_data("spatial", "coast.qs")) |> 
    sf::st_simplify(dTolerance = 100) |> 
    terra::vect() |> 
    terra::lines()
  terra::sbar(10000)
  terra::plot(mpa, xlim = terra::ext(mpa)[1:2], ylim = terra::ext(mpa)[3:4])
}

#### Get tagging locations 
# Get capture locations for relevant individuals 
tagsf <- 
  skateids |> 
  filter(individual_id %in% unique(iteration_patter$individual_id)) |> 
  select(lon = long_tag_capture, lat = lat_tag_capture) |>
  as.matrix() |> 
  terra::project(from = "EPSG:4326", to = terra::crs(mpa)) |>
  as.data.frame() |> 
  sf::st_as_sf(coords = c("V1", "V2"), crs = terra::crs(mpa))
# Write to file
write.csv(sf::st_coordinates(tagsf), 
          here_data_raw("movement", "tagging-sites.csv"), 
          row.names = FALSE)

#### Get receiver locations
moorings |>
  select(receiver_x, receiver_y) |> 
  as.data.frame() |> 
  write.csv(here_data_raw("movement", "moorings.csv"), 
            row.names = FALSE)

#### (deprecated) Plot study area via ggplot2
if (FALSE) {
  tic()
  png(here_fig("study-area.png"), 
      height = 5, width = 4, units = "in", res = 800)
  op <- options(terra.pal = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(256)), 0.75))
  # terra::plot(terra::rast(here_data("spatial", "bathy.tif")))
  bb   <- terra::ext(bathy) - 6.5e3
  xlim <- bb[1:2]
  ylim <- bb[3:4]
  p <- 
    ggplot_maps(data.table(mapfile = here_data("spatial", "bathy.tif"), row = 1, column = 1), 
                xlim = xlim,
                ylim = ylim,
                zlim = c(0, 350),
                mask = TRUE, 
                png_args = NULL) 
  p + 
    geom_sf(data = tagsf, shape = 8, size = 0.9, colour = "darkgreen") + 
    coord_sf(xlim = xlim, ylim = ylim) + 
    theme(panel.grid.major = element_line(colour = "#0000FF", linewidth = 0.025))
  options(op)
  dev.off()
  toc()
}


###########################
#### Example UDs from each algorithm (siam-linux20)

if (on_server() & FALSE) {
  
  #### Define mapfiles
  algorithms <- c("COA", "RSP", "AC", "DC", "ACDC")
  mapfiles <-
    data.table(algorithm = factor(algorithms, algorithms),
               mapindex  = seq_len(length(algorithms)),
               mapfile   = c("coa/ac/1/ud/spatstat/h/ud.tif", 
                             "rsp/ac/1/ud/dbbmm/ud.tif", 
                             "patter/ac/1/ud/spatstat/h/ud.tif", 
                             "patter/dc/1/ud/spatstat/h/ud.tif", 
                             "patter/acdc/1/ud/spatstat/h/ud.tif"))
  
  #### Collate data.table for mapping
  mapfiles <- 
    CJ(individual_id = c(25, 27, 35), # c(13, 20, 25, 27, 29, 35, 36, 38),
       mapindex      = mapfiles$mapindex
    ) |>  
    mutate(row    = factor(individual_id, levels = sort(unique(individual_id))), 
           column = mapfiles$algorithm[match(mapindex, mapfiles$mapindex)],
           mapfile = mapfiles$mapfile[match(mapindex, mapfiles$mapindex)],
           mapfile = file.path("data", "output", "analysis", individual_id, "04-2016", mapfile)) |> 
    select(row, column, mapfile) |> 
    as.data.table()
  
  #### File copy onto siam-linux20
  # data/output/analysis/25/04-2016/coa
  # data/output/analysis/25/04-2016/rsp
  # data/output/analysis/27/04-2016/coa
  # data/output/analysis/27/04-2016/rsp
  # data/output/analysis/35/04-2016/coa/
  # data/output/analysis/35/04-2016/rsp/
  # data/output/analysis/36/04-2016/coa/
  # data/output/analysis/36/04-2016/rsp/
  
  #### Make maps 
  ggplot_maps(mapfiles, 
              png_args = list(filename = here_fig("analysis", "map-examples.png"), 
                              height = 4, width = 4, units = "in", res = 800))
  
  pdf(here_fig("analysis", "map-examples.pdf"))
  p <- ggplot_maps(mapfiles)
  print(p)
  dev.off()
  
}


###########################
#### Overall ACDC UD (siam-linux20)

if (on_server()) {
  
  #### Estimate overall UD
  # List tif files
  mapfiles <- 
    iteration_patter |> 
    filter(dataset == "ACDC" & sensitivity == "Best") |> 
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
  
  #### Get angling records
  # Get all angling records (download from 28/11/2024)
  # (We expect ! NAs introduced by coercion warning here)
  cr <- 
    here_data_raw("movement", "skatespotter", "data (15).csv") |> 
    read.csv() |> 
    select(individual_id, date = date_captured, lon = longitude, lat = latitude) |> 
    mutate(date = as.Date(date), 
           lon = as.numeric(lon), 
           lat = as.numeric(lat)) |>
    filter(!is.na(lon) & !is.na(lat)) |> 
    filter(lon != 0 & lat != 0) |>
    as.data.table()
  # Check temporal distribution
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
  # Exclude locations on land 
  crxy <- crxy[!is.na(terra::extract(bathy_5m, cbind(crxy$x, crxy$y))$map_value), ]
  crsf <- 
    sf::st_as_sf(crxy, coords = c("x", "y"), crs = terra::crs(mpa))
  
  #### Quick plot of UD
  terra::plot(ud)
  points(crxy$x, crxy$y)
  hr <- map_hr_home(ud, .add = TRUE)
  poly <- terra::as.polygons(hr == 1)
  poly <- poly[poly[[1]] == 1]
  poly <- poly |> sf::st_as_sf()
  
  #### Plot UD with tagging & angling records
  pdf(here_fig("analysis", "map-overall.pdf"), 
      height = 3.5, width = 2.5)
  # Define plot (incl. limits)
  p <- 
    ggplot_maps(data.table(mapfile = ud.tif, row = 1, column = 1), 
                png_args = NULL)
  # Update plot & fix limits to those defined above
  p <- 
    p + 
    geom_sf(data = poly, fill = NA) + 
    geom_sf(data = crsf, shape = 21, size = 0.001, linewidth = 0, colour = "purple", alpha = 0.2) + 
    geom_sf(data = tagsf, shape = 8, size = 0.9, colour = "darkgreen") + 
    coord_sf(xlim = p$coordinates$limits$x, ylim = p$coordinates$limits$y) + 
    theme(panel.grid.major = element_line(colour = scales::alpha("dimgrey", 0.5), linewidth = 0.0001))
  print(p)
  dev.off()
  
}


###########################
#### QGIS inputs

if (FALSE) {
  
  #### Colour hexs
  # Land fill: #6969694C
  # Land border: blank (thin)
  scales::alpha("dimgrey", 0.3)
  # MPA open border: #0000FF80
  scales::alpha("blue", 0.5)
  # MPA closed border: #FF0000CC
  scales::alpha("#FF000080", 0.8)
  # Tagging locations: #006400FF
  scales::alpha("darkgreen", 1.0)
  
  #### Bathymetry colour bar
  breaks <- seq(-350, 0, by = 1); length(breaks)
  labels <- breaks
  labels[!(breaks %in% c(-300, -200, -100, 0))] <- ""
  png(here_fig("study-site-legend.png"),
      height = 5, width = 5, units = "in", res = 600)
  ggplot(data.frame(depth = 00:350)) + 
    geom_line(aes(depth, 1, colour = depth)) + 
    scale_color_gradientn(
      breaks = breaks, 
      labels = labels,
      colours = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(256)),
      limits = c(-350, 0),
      na.value = "grey", 
      guide = guide_colourbar(ticks = FALSE)) +
    theme(
      legend.text = element_text(margin = margin(t = 0, r = 3, b = 0, l = -1)) 
    )
  dev.off()
  
}


###########################
#### Algorithm sensitivity (siam-linux20)

if (on_server()) {
  
  id  <- 27        # 25, 27, 36
  mon <- "04-2016"
  
  #### Define mapfiles
  mapfiles <- rbindlist(
    list(
      # COA mapfiles
      # * Manually copy data/output/analysis/36/04-2016/coa onto server
      iteration_coa |> 
        filter(individual_id == id & month_id == mon) |> 
        mutate(algorithm = "COA", 
               sensitivity = factor(delta_t,
                                    c("2 days", "1 day", "3 days"), 
                                    labels = c("Best", "AC(-)", "AC(+)")), 
               mapfile = file.path(folder_ud, "spatstat", "h", "ud.tif")) |> 
        select(row = sensitivity, column = algorithm, mapfile) |>
        as.data.table(),
      
      # RSP mapfiles
      # * Manually copy data/output/analysis/36/04-2016/coa onto server
      iteration_rsp |> 
        filter(individual_id == id & month_id == mon)  |> 
        mutate(algorithm = "RSP", 
               sensitivity = factor(er.ad,
                                    c(500, 250, 750), 
                                    labels = c("Best", "AC(-)", "AC(+)")), 
               mapfile = file.path(folder_ud, "dbbmm", "ud.tif")) |> 
        select(row = sensitivity, column = algorithm, mapfile) |>
        as.data.table(),
      # Patter mapfiles
      iteration_patter |> 
        filter(individual_id == id & month_id == mon) |>
        mutate(mapfile = file.path(folder_ud, "spatstat", "h", "ud.tif")) |> 
        select(row = sensitivity, column = dataset, mapfile) |>
        as.data.table()
    )
  ) |>
    # (optional) Swap rows and columns
    mutate(row0 = row, column0 = column, 
           row = column0, column = row0) |> 
    as.data.table()
  
  #### Make maps
  # If sensitivity (row) by algorithm (column), use: height = 6, width = 3, 
  # If algorithm (row) by sensitivity (column), use: height = 5, width = 5
  head(mapfiles)
  ggplot_maps(mapdt = mapfiles, 
              png_args = 
                list(filename = here_fig(
                  "analysis", 
                  paste0("map-sensitivity-", id, "-", mon, ".png")
                ), 
                height = 5, width = 5, units = "in", res = 800))
  
}


###########################
###########################
#### Residency 

if (FALSE) {
  
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
  pdf(here_fig("analysis", "residency-best.pdf"), 
      height = 6 * 7.25/10, width = 7.25)
  lw <- 0.6 # residency-best
  # lw <- 1   # residency
  residency |>
    filter(sensitivity == "Best") |>
    # filter(zone == "total") |> 
    as_tibble() |> 
    ggplot() + 
    geom_point(aes(
      month, time, 
      colour = individual_id, 
      # shape = sensitivity, 
      alpha = if_else(sensitivity == "Best", 1, 0.5), 
      size = if_else(sensitivity == "Best", lw, 0.5)
    )) + 
    geom_line(aes(
      month, time, 
      colour = individual_id, 
      group = interaction(individual_id, sensitivity), 
      alpha = if_else(sensitivity == "Best", 1, 0.5), 
      size = if_else(sensitivity == "Best", lw, 0.5)
    )) +
    # geom_smooth(aes(month, time), lwd = 1.5, col = "black", se = TRUE) + 
    scale_alpha_identity() +
    scale_size_identity() +
    # scale_shape_manual(values = c(20, 17, 15, 3, 4, 8, 13)) + 
    # scale_x_date(labels = scales::date_format("%b-%y")) +
    # For residency-best, tweak labels
    scale_x_date(
      breaks = as.Date(c("2016-04-01", "2016-07-01", "2016-10-01", 
                         "2017-01-01", "2017-04-01")),
      labels = c("", "Jul-16", "", "Jan-17", "")
    ) + 
    scale_y_continuous(expand = c(0, 0.05), limits = c(0, 1)) + 
    geom_hline(data = res_null, aes(yintercept = time, colour = NULL), linetype = "dashed") +
    xlab("Time (months)") + 
    ylab("Residency") + 
    guides(
      colour = guide_legend(order = 1, title = "Individual"),
      # shape = guide_legend(order = 2, title = "Parameterisation"),
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
  
}


#### End of code. 
###########################
###########################