if (!patter:::os_linux() | Sys.getenv("JULIA_SESSION") == "FALSE") {

ggplot_maps <- function(mapdt,
                        xlim = NULL, 
                        ylim = NULL,
                        zlim = NULL,
                        mask = FALSE,
                        simplify = TRUE,
                        grid = terra::rast(here_data("spatial", "ud-grid.tif")),
                        mpa = qreadvect(here_data("spatial", "mpa.qs")),
                        coast = qs::qread(here_data("spatial", "coast.qs")), 
                        moorings = qs::qread(here_data("input", "mefs", "moorings.qs")),
                        png_args = NULL) {
  
  tictoc::tic()
  on.exit(tictoc::toc(), add = TRUE)
  
  # Validate inputs
  # mapdt must be a dataframe with the following columns:
  # * row, column : row/column identifiers
  # * mapfile     : the path to the UD
  stopifnot(c("row", "column", "mapfile") %in% colnames(mapdt))
  n_row <- length(unique(mapdt$row))
  n_col <- length(unique(mapdt$column))
  
  # (optional) Simplify coastline for improved speed 
  # * Full resolution:  160 s for 12 plots
  # * dTolerance = 500: 3 s for 12 plots (& reasonable approximation of coast)
  # * dTolerance = 100: 3 s for 12 plots (& good approximation of coast)
  if (simplify) {
    coast <- sf::st_simplify(coast, dTolerance = 100) 
  } else {
    coast <- sf::st_as_sf(coast)
  }

  # (optional) Prepare MPA
  mpa_poly <- sf::st_as_sf(mpa)
  
  # Define a blank data.frame
  # * This is used to handle blank panels (for which we couldn't compute the UD)
  blank <- as.data.frame(grid, xy = TRUE, na.rm = TRUE) |> as.data.table()
  colnames(blank) <- c("x", "y", "map_value")
  blank[, map_value := NA_real_]
  blank[, col := NA_character_]
  
  # Build a mapdata data.frame
  # * This includes, for each panel (row/column), the raster coordinates & values
  mapdata <- lapply(split(mapdt, seq_row(mapdt)), function(d) {
    if (file.exists(d$mapfile)) {
      r <- terra::rast(d$mapfile)
      if (mask) {
        r <- terra::mask(r, coast, inverse = TRUE, touches = FALSE)
      }
      # Get map coordinates via as.data.frame() or terra::spatSample() 
      rdt <- as.data.frame(r, xy = TRUE) 
      rdt <- rdt[!is.na(rdt[, 3]), ]
      setDT(rdt)
      colnames(rdt) <- c("x", "y", "map_value")
      # Define colours
      # * facet_wrap() forces the same zlim across all plots
      # * For individual colours, between min and max, we have to:
      # - Build the col column here
      # - Use scale_fill_identity()
      cols <- getOption("terra.pal")
      if (is.null(zlim)) {
        zlim <- range(rdt$map_value)
      }
      ints <- seq(zlim[1], zlim[2], length.out = length(cols))
      cols <- data.table(int = ints, 
                         col = cols)
      rdt[, col := cols$col[findInterval(rdt$map_value, cols$int)]]
      stopifnot(all(!is.na(rdt$col)))
    } else {
      rdt <- copy(blank)
    }
    # Link data.tables (e.g., row/column)
    # cbind(d, rdt)
    rdt[, c("row", "column") := .(d$row, d$column)]
  }) 
  if (length(mapdata) == 1L) {
    mapdata <- mapdata[[1]]
  } else {
    mapdata <- rbindlist(mapdata)
  }
  
  # Update coast with gg mapping
  # * We duplicate coast for each raster
  # * We use a more transparent colour for panels with convergence failures
  # * We use a blank panel otherwise (for NA factor levels)
  mapdata[, key := paste(row, column)]
  coast <- 
    lapply(unique(mapdata$key), function(key) {
      coast$key <- key
      coast
    }) |> 
    dplyr::bind_rows() |> 
    mutate(row = mapdata$row[match(key, mapdata$key)], 
           column = mapdata$column[match(key, mapdata$key)], 
           col = mapdata$col[match(key, mapdata$key)], 
           alp = ifelse(is.na(col), 0.075, 0.3),
           col = scales::alpha("dimgrey", alp),
           key = NULL)
  
  # Update mpa with gg mapping
  mpa_poly <- 
    lapply(unique(mapdata$key), function(key) {
      mpa_poly$key <- key
      mpa_poly
    }) |> 
    dplyr::bind_rows() |> 
    mutate(row = mapdata$row[match(key, mapdata$key)], 
           column = mapdata$column[match(key, mapdata$key)],
           blank = is.na(mapdata$col[match(key, mapdata$key)]), 
           # col = ifelse(blank, NA, col),
           alp = ifelse(blank, 0.25, 0.5),
           col = scales::alpha(col, alp),
           key = NULL)
  
  # Update moorings with gg mapping
  # * This is necessary so that we only add receivers to non-blank panels
  moorings <- 
    lapply(unique(mapdata$key), function(key) {
    m <- copy(moorings)
    m[, key := key]
  }) |> 
    dplyr::bind_rows() |> 
    mutate(row = mapdata$row[match(key, mapdata$key)], 
           column = mapdata$column[match(key, mapdata$key)]) |>
    filter(key %in% mapdata$key) |> 
    as.data.table()
  
  # Define map limits
  if (is.null(xlim)) {
    # xlim <- terra::ext(grid)[1:2]
    xlim <- terra::ext(mpa)[1:2]
  }
  if (is.null(ylim)) {
    # ylim <- terra::ext(grid)[3:4]
    ylim <- terra::ext(mpa)[3:4]
  }
  
  # Build ggplot
  p <- 
    mapdata |> 
    ggplot() + 
    theme_bw() + 
    geom_raster(aes(x = x, y = y, fill = col)) +
    scale_fill_identity() + 
    # scale_fill_gradientn(colours = getOption("terra.pal"), na.value = "white") +
    geom_point(data = moorings, aes(receiver_x, receiver_y), shape = 4, size = 0.2, stroke = 0.2) + 
    geom_sf(data = coast, aes(fill = I(col)), linewidth = 0.15) + 
    geom_sf(data = mpa_poly, aes(color = I(col)), fill = NA, linewidth = 0.25) + 
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) + 
    xlab("") + ylab("") 
  if (nrow(mapdt) > 1L) {
    p <- 
      p +  
      facet_grid(row ~ column, drop = FALSE) + 
      theme(axis.text = element_blank(), 
            axis.ticks = element_blank(), 
            panel.grid = element_blank(), 
            # strip.text = element_blank(), 
            legend.position = "none")
  }
  
  # Return ggplot 
  if (!is.null(png_args)) {
    do.call(grDevices::png, png_args)
    print(p)
    dev.off()
    return(invisible(NULL))
  } else {
    return(p)
  }

}

}
