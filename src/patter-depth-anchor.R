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
                  main = paste0(sort(d$receiver_id), collapse = ", "))
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
    # id <- ids[1]
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
    
    # Define receiver_id_next_key (used for matching)
    obs[, receiver_id_next_key := 
          lapply(obs$receiver_id_next, acs_setup_receiver_key) |> unlist()]
    obs$receiver_id_next_key[obs$receiver_id_next_key == ""] <- NA_character_
    
    # Identify distinct container/depth combinations for selected individual
    obs |> 
      select("receiver_id_next_key", "depth") |> 
      distinct() |> 
      filter(!is.na(receiver_id_next_key)) |>
      as.data.table()
    
  }) |> 
    plyr::compact() |>
    rbindlist() |>
    distinct() |> 
    as.data.table()
  
}

#' @title AC* set up: detection container cells
#' @description
#' For each unique detection container (receiver(s)--depth combination), this function identifies the coordinates of valid cells. 

acs_setup_container_cells <- function(.dlist, .containers, .rcd) {
  
  # Check user inputs 
  check_dlist(.dlist = .dlist, 
              .spatial = c("bathy", "bset"))
  bathy <- .dlist$spatial$bathy
  bset  <- .dlist$spatial$bset
  
  # Loop over unique detection containers & build valid cell lists
  # * We save lists to file in ./data/input/containers/{container}/{depth}.parquet
  pbapply::pblapply(split(.rcd, .rcd$receiver_id_next_key), function(d) {
    
    # d <- split(.rcd, .rcd$receiver_id_next_key)[[1]]
    message(.rcd$index[1])
    
    # Identify cells within the container
    container <- .containers[[d$receiver_id_next_key[1]]]
    bathy_in_container <- terra::crop(bathy, container)
    bathy_in_container <- terra::mask(bathy_in_container, container)
    cells_in_container <- 
      bathy_in_container |>
      terra::as.data.frame(cell = TRUE, xy = TRUE) |> 
      mutate(digi = terra::extract(bset, cell)) |> 
      as.data.table()
    
    # Loop over depth observations & identify valid cells for each depth 
    cells_in_container_by_depth <- 
      pbapply::pblapply(split(d, collapse::seq_row(d)), function(.d) {
        
        # Define valid cells
        # .d <- d[1, ]
        dt <-
          cells_in_container |> 
          copy() |>
          # Define spatially explicit shallow and deep depth limits 
          calc_depth_envelope(.depth = .d$depth, 
                              .calc_depth_error = calc_depth_error) |>
          # Identify valid cells within container
          filter(depth > shallow_2 & depth < deep) |>
          select("x", "y") |>
          as.data.table()
        
        # Write to file 
        stopifnot(nrow(dt) > 0L) 
        arrow::write_parquet(dt, 
                             here_data("input", "containers", 
                                       .d$receiver_id_next_key,
                                       paste0(.d$depth, ".parquet")))
        
      })
  })
  invisible(NULL)
}
