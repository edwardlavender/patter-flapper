ggplot_maps <- function(mapdt,
                        xlim = NULL, 
                        ylim = NULL, 
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
  
  # Build a mapdata data.frame
  # * This includes, for each panel (row/column), the raster coordinates & values
  mapdata <- lapply(split(mapdt, seq_row(mapdt)), function(d) {
    r <- terra::rast(d$mapfile)
    # Get map coordinates via as.data.frame() or terra::spatSample() 
    rdt <- as.data.frame(r, xy = TRUE, na.rm = TRUE) |> as.data.table()
    colnames(rdt) <- c("x", "y", "map_value")
    cbind(d, rdt)
  }) |> rbindlist()

  # Build ggplot
  p <- 
    mapdata |> 
    ggplot() + 
    theme_bw() + 
    geom_raster(aes(x = x, y = y, fill = map_value)) +
    scale_fill_gradientn(colours = getOption("terra.pal"), na.value = "white") +
    geom_sf(data = coast, fill = scales::alpha("dimgrey", 0.5)) + 
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) + 
    xlab("") + ylab("") + 
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          panel.grid = element_blank(), 
          strip.text = element_blank(), 
          legend.position = "none") + 
    facet_wrap(~row + column, nrow = n_row, ncol = n_col)
  
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
