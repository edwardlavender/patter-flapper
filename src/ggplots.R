if (!patter:::os_linux()) {

ggplot_maps <- function(mapdt,
                        xlim = NULL, 
                        ylim = NULL,
                        grid = terra::rast(here_data("spatial", "ud-grid.tif")),
                        coast = qs::qread(here_data("spatial", "coast.qs")), 
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
  coast <- sf::st_simplify(coast, dTolerance = 100) |> sf::st_geometry()
  
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
      # Get map coordinates via as.data.frame() or terra::spatSample() 
      rdt <- as.data.frame(r, xy = TRUE, na.rm = TRUE) |> as.data.table()
      colnames(rdt) <- c("x", "y", "map_value")
      # Define colours
      # * facet_wrap() forces the same zlim across all plots
      # * For individual colours, between min and max, we have to:
      # - Build the col column here
      # - Use scale_fill_identity()
      cols <- getOption("terra.pal")
      ints  <- seq(min(rdt$map_value), max(rdt$map_value), length.out = length(cols))
      cols <- data.table(int = ints, 
                         col = cols)
      rdt[, col := cols$col[findInterval(rdt$map_value, cols$int)]]
      stopifnot(all(!is.na(rdt$col)))
    } else {
      rdt <- copy(blank)
    }
    # Link data.tables (e.g., row/column)
    cbind(d, rdt)
  }) |> rbindlist()
  
  # Define map limits
  if (is.null(xlim)) {
    xlim <- terra::ext(grid)[1:2]
  }
  if (is.null(ylim)) {
    ylim <- terra::ext(grid)[3:4]
  }

  # Build ggplot
  p <- 
    mapdata |> 
    ggplot() + 
    theme_bw() + 
    geom_raster(aes(x = x, y = y, fill = col)) +
    scale_fill_identity() + 
    # scale_fill_gradientn(colours = getOption("terra.pal"), na.value = "white") +
    geom_sf(data = coast, fill = scales::alpha("dimgrey", 0.5)) + 
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) + 
    xlab("") + ylab("") + 
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          panel.grid = element_blank(), 
          # strip.text = element_blank(), 
          legend.position = "none") + 
    facet_wrap(~row + column, drop = TRUE, nrow = n_row, ncol = n_col)
  
  # Return ggplot 
  if (!is.null(png)) {
    do.call(grDevices::png, png_args)
    print(p)
    dev.off()
    return(invisible(NULL))
  } else {
    return(p)
  }

}

}
