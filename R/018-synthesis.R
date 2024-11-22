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
acoustics_raw     <- qs::qread(here_data("input", "mefs", "acoustics.qs"))
archival_raw      <- qs::qread(here_data("input", "mefs", "archival.qs"))
acoustics_by_unit <- qs::qread(here_data("input", "acoustics_by_unit.qs"))
archival_by_unit  <- qs::qread(here_data("input", "archival_by_unit.qs"))
behaviour_by_unit <- qs::qread(here_data("input", "behaviour_by_unit.qs"))
recaps            <- readRDS(here_data_raw("movement", "recaptures_processed.rds"))


###########################
###########################
#### Visualise time series 

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
    geom_vline(data = recaps, aes(xintercept = timestamp), 
               colour = "red", alpha = 0.5,  linewidth = 0.35) + 
    geom_text(data = acoustics,
              aes(timestamp, 2.5, label = "|"),
              colour = "black", size = 2, vjust = 0.5) + 
    scale_x_datetime(labels = scales::date_format("%d"), 
                     breaks = xticks, 
                     expand = c(0, 0)) +
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
  
}


###########################
###########################
#### Convergence



###########################
###########################
#### Patterns of space use



###########################
###########################
#### Residency 



#### End of code. 
###########################
###########################
