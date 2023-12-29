#' @title Internal patter functions

utils.add::load_internal_functions("patter")

#' @title AC* set up: define receiver keys
#' @description This function strings together receiver combinations. 

acs_setup_receiver_key <- function(.receiver) {
  paste(sort(unlist(.receiver)), collapse = "-")
}

#' @title AC* set up: detection container receiver combinations
#' @title This function defines the combinations of receivers that recorded detections simultaneously. 

acs_setup_containers_rc <- function(.dlist, .split = NULL) {
  check_dlist(.dlist = .dlist, 
              .dataset = "acoustics")
  acoustics <- copy(.dlist$data$acoustics)
  check_names(acoustics, c(.split, "timestamp", "receiver_id"))
  if (is.null(.split)) {
    if (rlang::has_name(acoustics, "individual_id") & 
        length(unique(acoustics$individual_id)) > 1L) {
      warn("The acoustic data contains an `individual_id` column with multiple individuals but `.split = NULL`.")
    }
    acoustics[, individual_id := 1L]
  } else {
    acoustics[, individual_id := acoustics[[.split]]]
  }
  acoustics |>
    # Round time steps 
    mutate(timestamp = lubridate::round_date(.data$timestamp, "2 mins")) |>
    group_by(.data$individual_id, .data$receiver_id, .data$timestamp) |>
    slice(1L) |>
    ungroup() |>
    # List receivers with detections at each time step
    group_by(.data$individual_id, timestamp) |>
    arrange(.data$receiver_id, .by_group = TRUE) |>
    summarise(receiver_id = list(unique(receiver_id)), 
              receiver_key = acs_setup_receiver_key(.data$receiver_id)) |> 
    ungroup() |>
    # Identify unique receiver combinations (across all individuals & time steps)
    filter(!duplicated(.data$receiver_key)) |> 
    mutate(index = row_number()) |> 
    select("index", "receiver_id", "receiver_key") |> 
    as.data.table() 
}

#' @title AC* setup up: detection containers
#' @description This function builds a named `list` of detection containers, one for each unique combination of receivers that recorded detections simultaneously (`.rc`).

acs_setup_detection_containers <- function(.dlist, .rc, .plot = FALSE) {
  
  # Check user inputs 
  check_dlist(.dlist = .dlist, 
              .dataset = "moorings", 
              .spatial = "bathy")
  moorings <- copy(dlist$data$moorings)
  check_names(moorings, c("receiver_id",
                          "receiver_x", "receiver_y", 
                          "receiver_range"))
  check_names(.rc, "receiver_id")
  
  # Identify relevant data in dlist
  gamma <- moorings$receiver_range[1]
  moorings <- 
    moorings |> 
    select("receiver_id", "receiver_x", "receiver_y") |>
    as.data.table()
  crs <- terra::crs(.dlist$spatial$bathy)
  
  # Define a list of detection containers
  containers <- pbapply::pblapply(split(.rc, collapse::seq_row(.rc)), function(d) {
    
    # Define receiver coordinates
    # d <- .rc[361, ]
    d <- 
      data.table(receiver_id = unlist(d$receiver_id)) |>
      left_join(moorings |> 
                  filter(receiver_id %in% receiver_id), 
                by = "receiver_id") |> 
      as.data.table()
    
    # Define containers
    containers <- 
      lapply(collapse::seq_row(d), function(i) {
        cbind(d$receiver_x[i], d$receiver_y[i]) |> 
          terra::vect(crs = crs) |> 
          terra::buffer(width = gamma)
      })
    container <- spatIntersect(containers)
    
    # Optionally plot
    if (.plot) {
      cols <- seq_len(length(containers))
      containers <- do.call(rbind, containers)
      terra::plot(containers, 
                  border = cols, 
                  main = acs_setup_receiver_key(d$receiver_id))
      text(d$receiver_x, d$receiver_y, d$receiver_id, 
           font = 2, cex = 0.5, col = cols)
      terra::plot(container, add = TRUE, col = scales::alpha("lightgrey", 0.5))
      readline("Press [Enter] to continue...")
    }
    
    # Return outputs
    container
    
})
  
  names(containers) <- .rc$receiver_key
  containers
  
}

#' @title AC* set up: detection container receiver/depth combinations
#' @description This function defines the unique combinations of receivers and depth observations at the moment of detection.

acs_setup_containers_rcd <- function(.dlist) {
  
  # * Pair time series for each individual & join
  # * Include container IDs
  # * Identify unique receiver/depth combinations
  
  # Check user inputs
  check_dlist(.dlist = dlist, 
              .dataset = c("acoustics", "archival"))
  acoustics <- .dlist$data$acoustics
  archival  <- .dlist$data$archival
  check_names(acoustics, "individual_id")
  check_names(archival, "individual_id")
  
  # Define individuals with acoustic and archival data 
  acc_ids <- unique(acoustics$individual_id)
  ids     <- acc_ids[acc_ids %in% unique(archival$individual_id)]
  
  # Define (unique) container/depth combinations (across all ids)
  pbapply::pblapply(ids, function(id) {
    
    # Isolate individual-specific acoustic & archival data 
    # (id <- ids[2])
    acc   <- acoustics[individual_id == id, ]
    arc   <- archival[individual_id == id, ]
    dlist <- pat_setup_data(.acoustics = acc, 
                            .archival = arc)
    
    # Collate observations with .trim = TRUE
    # * Use tryCatch() to handle individuals without overlapping time series 
    obs <- tryCatch(
      pf_setup_obs(.dlist = dlist, 
                   .trim = TRUE, 
                   .step = "2 mins", 
                   .mobility = 500,
                   .receiver_range = 750), 
      error = function(e) e)
    if (inherits(obs, "error")) {
      warn(obs$message)
      return(NULL)
    }
    
    # Focus on detections
    # * This is essential to identify the correct combinations
    # * ... of receiver_id_next_key & depths align 
    obs <- obs[detection == 1L, ]
    # Use next depths to align with receiver_id_next_key
    obs[, depth := lead(depth)]
    
    # Define receiver_id_next_key (used for matching)
    obs[, receiver_id_next_key := 
          lapply(obs$receiver_id_next, acs_setup_receiver_key) |> unlist()]
    obs$receiver_id_next_key[obs$receiver_id_next_key == ""] <- NA_character_
    
    # Identify distinct container/depth combinations for selected individual
    # * Retain individual_id and timestamp to facilitate debugging 
    obs[, individual_id := id]
    obs |> 
      select("individual_id", "timestamp", "receiver_id_next_key", "depth") |> 
      distinct(receiver_id_next_key, depth, .keep_all = TRUE) |> 
      filter(!is.na(receiver_id_next_key)) |>
      as.data.table()
    
  }) |> 
    plyr::compact() |>
    rbindlist() |>
    distinct(receiver_id_next_key, depth, .keep_all = TRUE) |> 
    arrange(receiver_id_next_key, depth) |> 
    mutate(index = rleid(receiver_id_next_key)) |> 
    select("index", "individual_id", "timestamp", "receiver_id_next_key", "depth") |> 
    as.data.table()
  
}

#' @title AC* set up: detection container cells
#' @description For each unique detection container (receiver(s)--depth combination), this function identifies the coordinates of valid cells. 

acs_setup_container_cells <- function(.dlist, .containers, .rcd, 
                                      .cl = NULL) {
  
  # Check user inputs 
  check_dlist(.dlist = .dlist, 
              .spatial = c("bathy", "bset"))
  bathy <- .dlist$spatial$bathy
  bset  <- .dlist$spatial$bset
  
  # Loop over unique detection containers & build valid cell lists
  pbapply::pblapply(split(.rcd, .rcd$receiver_id_next_key), function(d) {
    
    # d <- split(.rcd, .rcd$receiver_id_next_key)[[19]]
    message(d$index[1])
    
    # Define outfiles
    # * We save lists to file in ./data/input/containers/{container}/{depth}.parquet
    d[, outfile := here_data("input", "containers", 
                             receiver_id_next_key,
                             paste0(depth, ".parquet"))]
    d[, outexists := file.exists(outfile)]
    d <- d[outexists == FALSE, ]
    if (nrow(d) == 0L) {
      return(NULL)
    }
    
    # Identify cells within the container
    container <- .containers[[d$receiver_id_next_key[1]]]
    bathy_in_container <- terra::crop(bathy, container)
    bathy_in_container <- terra::mask(bathy_in_container, container)
    if (FALSE) {
      # Plot container 
      terra::plot(container)
      terra::plot(bathy_in_container, add = TRUE)
      terra::lines(container)
      # Add coordinates for a selected receiver
      moorings |> 
        filter(receiver_id == 27) |> 
        select(receiver_easting, receiver_northing) |>
        points()
      # Check range in depths in container
      terra::global(bathy_in_container, "range", na.rm = TRUE)
    }
    cells_in_container <- 
      bathy_in_container |>
      terra::as.data.frame(xy = TRUE) |> 
      mutate(digi = terra::extract(bset, cbind(x, y))) |> 
      as.data.table()
    
    # Loop over depth observations & identify valid cells for each depth 
    # * This code is implemented in parallel
    # * We implement this inner loop in parallel b/c it does not require SpatRasters
    # * (which cannot be serialised over the connection without wrapping/unwrapping)
    cells_in_container_by_depth <- 
      cl_lapply(.x = split(d, collapse::seq_row(d)), 
                .cl = .cl, .use_chunks = TRUE, 
                .fun = function(.d) {
        
        # Define valid cells
        # .d <- d[1, ]
        print(.d)
        dt <- 
          cells_in_container |> 
          copy() |>
          # Define spatially explicit shallow and deep depth limits 
          mutate(cell_now = terra::cellFromXY(.dlist$spatial$bathy, cbind(x, y))) |>
          as.data.table() |>
          calc_depth_envelope(.obs = NULL, .t = 1L, .dlist = .dlist) |>
          # Identify valid cells within container
          filter(depth >= shallower & depth <= deep) |>
          select("x", "y") |>
          as.data.table()
        
        # Write to file 
        stopifnot(nrow(dt) > 0L) 
        arrow::write_parquet(dt, .d$outfile)
        NULL
        
      })
  })
  invisible(NULL)
}

#' @title AC* helper: ACDC container likelihood
#' @description This function filters particle proposals that are incompatible with container dynamics. When there is a long time before the next detection, we use acs_filter_container(). As we approach the next detection, we refine the representation of container dynamics with acs_filter_container_acdc(). This accounts for the observed depth at the time of the next acoustic detection. 

acs_filter_container_acdc <- function(.particles, .obs, .t, .dlist) {
  
  # Check user inputs
  if (.t == 1L) {
    check_dlist(.dlist = .dlist, 
                .algorithm = c("pos_detections"))
    stopifnot(length(unique(.obs$mobility)) == 1L)
  }
  
  # * Identify the time step of the next detection (if applicable)
  # * Identify the number of time steps before the next detection
  # * If there are more than N time steps before the next detection, 
  # * ... or there are very large numbers of location combinations
  # * ... we use the usual acs_filter_container() function only 
  # * ... (for speed & to avoid memory issues)
  # * If there are less than N time steps before the next detection, 
  # * ... we use this account for depth observations
  # * This reduces the burden of this filter, while accounting for the placement
  # * ... of valid locations by depth when it is most important 
  # * (... when we approach receiver(s))
  
  if (.t > 1 && .t < max(.obs$timestep)) {
    
    #### Define whether or not to implement the revised ACDC filter  
    
    # Start with do_acdc = TRUE
    do_acdc <- TRUE
    
    # Only implement the filter if the time gap until the next detection is short
    pos_detections <- .dlist$algorithm$pos_detections
    pos_detection  <- (pos_detections[pos_detections > .t])[1L]
    timegap <- pos_detection - .t
    if (timegap > 25L) {
      do_acdc <- FALSE
    }
    
    # Only implement the filter for valid depth observations
    if (do_acdc) {
      # Get the ID of the next receiver(s) and the corresponding depth observation
      container <- .obs$receiver_id_next_key[.t]
      depth     <- .obs$depth[pos_detection]
      if (is.na(depth)) {
        do_acdc <- FALSE
      }
    }
    
    # Only implement the filter if possible within vector memory
    if (do_acdc) {
      # Read valid locations at the next time step
      locs <- arrow::read_parquet(
        file.path("data", "input", "containers", 
                  container, 
                  paste0(depth, ".parquet")))
      # Calculate the number of calculations required
      nc <- nrow(.particles) * nrow(locs)
      # Handle integer overflow
      if (is.na(nc)) {
        nc <- Inf
      }
      # Only implement the filter if we can implement the calculations in memory (relatively quickly)
      if (nc > 80e6L) {
        do_acdc <- FALSE
      }
    }

    #### Implement the usual AC* filter
    if (!do_acdc) {
      .particles <- acs_filter_container(.particles = .particles, 
                                         .obs = .obs, 
                                         .t = .t, 
                                         .dlist = .dlist)
      
    #### Implement the revised ACDC filter that accounts for depth observations
    } else {
      # Calculate distances between particle samples & valid locations at the next detection
      # * Rows are particles
      # * Columns are future locations
      dist <- terra::distance(.particles |>
                                select("x_now", "y_now") |>
                                as.matrix(),
                              locs |>
                                as.matrix(),
                              lonlat = .dlist$par$lonlat)                    
      # Eliminates particles using distance threshold
      # * We eliminate particles that are not within a 
      # * ... reachable distance of at least one valid location
      .particles <- .particles[Rfast::rowsums(dist <= (.obs$mobility[.t] * timegap)) > 0L, ]
    }
  }
  
  #### Return outputs
  .particles
}