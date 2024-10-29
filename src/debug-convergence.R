#' Animate particles from the filter
#' @param .sim The data.table row
#' @param .map The SpatRaster
#' @param .moorings The moorings data.table
#' @param .start,.stop Integers that define the time steps of interest
#' @param .input The named list of arguments passed to pf_filter()
#' @param .output The named list of outputs from pf_filter()

ani <- function(.sim, .map, .moorings, .start = 1L, .end, .input, .output) {
  
  .moorings <- copy(.moorings)
  
  # Create directories
  frames <- here_fig("debug", .sim$index, "frames")
  dir.create(frames)
  mp4 <- here_fig("debug", .sim$index)
  dir.create(mp4)
  
  # Make frames
  pf_plot_xy(.map = .map, 
             .coord = .output$states, 
             .steps = .start:.end,
             .png = list(filename = frames),
             # .cl = 10L,
             .add_points = list(pch = ".", col = "red"),
             .add_layer = function(t) {
               
               # Add receivers
               text(.moorings$receiver_x, .moorings$receiver_y, .moorings$receiver_id, cex = 0.5)
               
               # Add detection containers
               .moorings[, receiver_gamma := 1750]
               cbind(.moorings$receiver_x, .moorings$receiver_y) |>
                 terra::vect() |>
                 terra::buffer(width = .moorings$receiver_gamma) |> 
                 terra::lines()
               
               # Add acoustic containers
               containers <- .input$.yobs$ModelObsAcousticContainer
               if (!is.null(containers)) {
                 cinfo <- containers[timestamp == .input$.timeline[t], ]
                 if (nrow(cinfo) > 0L) {
                   
                   # Colour receiver(s) with next detection
                   text(cinfo$receiver_x, cinfo$receiver_y, cinfo$sensor_id, cex = 0.5, col = "blue", font = 2)
                   
                   # Add time-specific acoustic containers
                   cbind(cinfo$receiver_x, cinfo$receiver_y) |>
                     terra::vect() |>
                     terra::buffer(width = cinfo$radius) |> 
                     terra::lines(col = "royalblue") 
                 }
               }
             })
  
  # Make animation 
  input   <- file_list(frames)
  output  <- file.path(mp4, "ani.mp4")
  av::av_encode_video(input, output, framerate = 5)
  
  # Open animation (on MacOS)
  system(paste("open", shQuote(output)))
  invisible(NULL)
}

#' Interactively debug convergence failures in the ACDC algorithm 

debug_acdc <- function(.map, .moorings, .step, .input, .output) {
  
  .moorings <- copy(.moorings)
  
  # Define time step of interest (i.e., a time step before a convergence failure)
  timeline <- .input$.timeline
  time     <- timeline[.step]
  
  # Examine observations 
  .yobs <- .input$.yobs
  time |> print()
  .yobs$ModelObsAcousticLogisTrunc[timestamp >= time, ] |> head() |> print()
  .yobs$ModelObsAcousticContainer[timestamp >= time, ] |> head() |> print()
  
  # Define map region for visualisation 
  containers <- .yobs$ModelObsAcousticContainer
  rxy        <- containers[timestamp == time, ]
  if (nrow(rxy) > 0L) {
    map_region <- terra::crop(.map, 
                              cbind(rxy$receiver_x, rxy$receiver_y) |>
                                terra::vect() |>
                                terra::buffer(width = max(rxy$radius)) |>
                                terra::ext())
  } else {
    map_region <- .map
  }
  
  # map_region <- map
  
  # Plot particles
  pf_plot_xy(.map = map_region, 
             .coord = .output$states, 
             .steps = .step:length(timeline), 
             .add_points = list(pch = ".", col = "red"),
             .add_layer = function(t) {
               
               # Add receivers
               text(.moorings$receiver_x, .moorings$receiver_y, .moorings$receiver_id, cex = 0.5)
               
               # Add detection containers
               .moorings[, receiver_gamma := 1750]
               cbind(.moorings$receiver_x, .moorings$receiver_y) |>
                 terra::vect() |>
                 terra::buffer(width = .moorings$receiver_gamma) |> 
                 terra::lines()
               
               # Colour receiver(s) with next detection
               cinfo <- containers[timestamp == timeline[t], ]
               text(cinfo$receiver_x, cinfo$receiver_y, cinfo$sensor_id, cex = 0.5, col = "blue", font = 2)
               
               # Add time-specific acoustic containers
               cbind(cinfo$receiver_x, cinfo$receiver_y) |>
                 terra::vect() |>
                 terra::buffer(width = cinfo$radius) |> 
                 terra::lines(col = "royalblue")
               
             }, 
             .prompt = TRUE)
  
  invisible(NULL)
  
}