#' Compute distances
dist_2d <- function(xy0, xy1) {
  sqrt((xy0[, 1] - xy1[, 1])^2 + (xy0[, 2] - xy1[, 2])^2)
}

#' Assemble (re)capture containers
# * .xinit is a named list of starting coordinates (forward, backward). Each element ia a data.table with one or more rows.
# * .radius is the min distance from the starting coordinates we must reach.

# TO DO
# * Add optional .map argument & define threshold from map_maxdim()

assemble_capture_containers <- function(.timeline, 
                                        .xinit = list(), 
                                        .radius,
                                        .mobility,
                                        .threshold = NULL) {
  
  # Check user inputs
  check_timeline(.timeline)
  check_named_list(.xinit)
  check_names(.xinit, c("forward", "backward"))
  directions <- c("forward", "backward")
  
  # Assemble containers (forward & backward)
  containers <- lapply(directions, function(.direction) {
    .assemble_capture_containers(.timeline  = .timeline, 
                                 .xinit     = .xinit, 
                                 .radius    = .radius,
                                 .mobility  = .mobility,
                                 .threshold = .threshold, 
                                 .direction = .direction)
  })
  
  # Return containers
  names(containers) <- directions
  containers
  
}

.assemble_capture_containers <- function(.timeline,
                                         .xinit = list(), 
                                         .radius,
                                         .mobility,
                                         .threshold = NULL,
                                         .direction = c("forward", "backward")) {
  
  # Define direction
  .direction <- match.arg(.direction)
  
  # Define container information for the 'ending' coordinate
  # * If direction = "forward", .xinit$backward defines the container coordinates
  # * If direction = "backward", .xinit$forward defines the container coordinates
  if (.direction == "forward") {
    cinfo <- .xinit$backward
  } else if (.direction == "backward") {
    cinfo <- .xinit$forward 
  }
  
  # Define container coordinates 
  # * cx, cy are container coordinates
  # * radius is container radius 
  if (nrow(cinfo) == 1L) {
    cx     <- cinfo$x
    cy     <- cinfo$y
    radius <- .radius
  } else {
    # Use the central point if multiple points provided
    # radius is defined as the maximum distance from the centre to a point 
    # (plus the .radius buffer)
    cx     <- mean(cinfo$x)
    cy     <- mean(cinfo$y)
    radius <- max(dist_2d(cbind(cx, cy), cbind(cinfo$x, cinfo$y)))
    radius <- radius + .radius
  }
  
  # Define sequence of radii
  # For .direction = "forward", radii shrink _forward_ in time
  # For .direction = "backward", radii shink _backward_ in time 
  # radii[1] or radii[T] = radius 
  radii <- radius + (1:length(.timeline) - 1) * .mobility
  if (.direction == "forward") {
    radii <- rev(radii)
  }
  
  # Build container 
  recap_container <- data.table(timestamp  = .timeline, 
                                obs        = 1L,
                                sensor_id  = 1L, 
                                capture_x  = cx,
                                capture_y  = cy,
                                radius     = radii)
  
  # Filter by .threshold 
  if (!is.null(.threshold)) {
    radius          <- NULL
    recap_container <- recap_container[radius <= .threshold, ]
  }
  
  recap_container
}

# Examples
# timeline <- seq(as.POSIXct("2016-01-01 00:00:00", tz = "UTC"), 
#                 as.POSIXct("2016-01-01 00:10:00", tz = "UTC"), 
#                 by = "2 mins")
# xinit    <- list(forward = data.table(x = 0, y = 0), 
#                  backward = data.table(x = 10, y = 10))
# assemble_capture_containers(.timeline = timeline, 
#                             .xinit = xinit, 
#                             .radius = 100,
#                             .mobility = 500,
#                             .threshold = 5000)

# Julia source wrapper
set_ModelObsCaptureContainer <- function() {
  JuliaCall::julia_source(here::here("Julia", "src", "ModelObsCaptureContainer.jl"))
}
