#' @title 2d Euclidean distances 

dist_2d <- function(x0, x1, y0, y1) {
  sqrt((x0 - x1)^2 + (y0 - y1)^2)
}

#' @title Internal {patter} helpers
cat_helper <- patter:::cat_init
call_start <- patter:::call_start
call_end   <- patter:::call_end
pb_init    <- patter:::pb_init
pb_tick    <- patter:::pb_tick
pb_close   <- patter:::pb_close

#' @title PF: run the backward pass
#' @description This is an implementation of the backward sampling algorithm.
#' * At each time step, we identify all combinations of cells and calculate densities;
#' * We assume that the likelihood evaluation is cheap, so it is quicker to calculate densities for all combinations of particles (using fast vectorised R code) than identify the unique cell of cell combinations, do the calculations for the subset of unique cells, and then match back on to a larger dataframe for resampling;
#' * This approach trades parallelisation over particles for vectorisation and is faster for moderate numbers of cores;
#' * This generates outputs like the forward run;
#' * We can use this with the usual functions;
#' * But we don't automatically get the trajectories;

pf_backward_sampler_2 <- function(.history, 
                                  .write_opts, 
                                  .progress = TRUE, .verbose = TRUE, .txt = "") {
  
  #### Check user inputs
  t_onset <- Sys.time()
  
  #### Set up messages
  cat_log <- cat_init(.verbose = .verbose)
  cat_log(call_start(.fun = "pf_forward", .start = t_onset))
  on.exit(cat_log(call_end(.fun = "pf_forward", .start = t_onset, .end = Sys.time())), add = TRUE)
  
  #### Set up loop
  # Number of time steps 
  n_step <- length(.history)
  # Read history 
  .history[[n_step]] <- arrow::read_parquet(.history[[n_step]])
  # Particle index 
  index <- collapse::seq_row(.history[[n_step]])
  # Global vars
  dens <- NULL
  # Progress bar
  pb <- pb_init(.n = n_step - 1L, .init = 0L, .progress = .progress)
  
  #### Run backward sampler
  for (t in n_step:2L) {
    
    cat_to_cf(paste0("... Time step ", t, ":"))
    pb_tick(.pb = pb, .t = n_step - (t - 1), .progress = .progress)
    
    #### Collect .history[[t]] and .history[[t - 1L]]
    .history[[t]][, index := index]
    .history[[t - 1L]] <- arrow::read_parquet(.history[[t - 1L]])
    .history[[t - 1L]][, index := index]
    
    #### Define .history[[t]]
    # For each particle at time t, we sample a previous particle
    h <-
      CJ(
        now_index = index, # .history[[t]]$index,
        past_index = index # .history[[t - 1L]]$index
      ) |> 
      lazy_dt() |> 
      left_join(.history[[t]] |> 
                  select("index", "cell_now", "x_now", "y_now"), 
                by = c("now_index" = "index")) |>
      left_join(.history[[t - 1]] |> 
                  select(index, 
                         cell_past = "cell_now", 
                         x_past = "x_now", 
                         y_past = "y_now"), 
                by = c("past_index" = "index")) |>
      mutate(dist = dist_2d(x0 = .data$x_past, x1 = .data$x_now, 
                            y0 = .data$y_past, y1 = .data$y_now),
             dens = dtruncgamma2(.data$dist)) |>
      group_by(.data$now_index) |>
      # slice_sample() doesn't work with .data$dens pronoun
      slice_sample(n = 1L, weight_by = dens) |>
      ungroup() |>
      mutate(timestep = t - 1L) |>
      select("timestep", 
             "cell_past", "x_past", "y_past", 
             "cell_now", "x_now", "y_now", 
             "dist", "dens") |>
      as.data.table()
    
    #### Record outputs
    arrow::write_parquet(h, file.path(.write_opts$sink, paste0(t, ".parquet")))
    .history[[t]] <- NA
    
    #### Move on
    .history[[t - 1L]] <- 
      h |> 
      select(cell_now = "cell_past", x_now = "x_past", y_now = "y_past") |>
      as.data.table()
    
  }
  pb_close(.pb = pb, .progress = .progress)
  
  #### Record outputs for .history[[t]]
  .history[[1L]] |> 
    mutate(timestep = 1L, 
           cell_past = NA_integer_, 
           x_past = NA_real_, 
           y_past = NA_real_, 
           dist = NA_real_, 
           dens = NA_real_) |> 
    select("timestep", 
           "cell_past", "x_past", "y_past",
           "cell_now", "x_now", "y_now", 
           "dist", "dens") |> 
    arrow::write_parquet(file.path(.write_opts$sink, paste0(1, ".parquet")))
  
  #### Return outputs 
  NULL
  
}

