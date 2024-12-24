#' Cat helpers

cat_line <- function() {
  cat("\n\n\n---------------------------------------------------------------\n")
}

cat_iter <- function(t) {
  msg("\n On row {t}...", .envir = environment())
}

cat_init <- function(t) {
  cat_line()
  cat_iter(t)
}

#' Create directories 
dirs.create <- function(paths) {
  sapply(paths, function(path) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
    }
  }) |> 
    suppressWarnings()
  invisible(NULL)
}

#' Copy directories (and contents)

dirs.copy <- function(from, to, cl) {
  
  # Assemble data.table and handle duplicates
  folders <- 
    data.table(from = from, to = to) |> 
    group_by(from) |> 
    slice(1L) |> 
    as.data.table()

  # Validate `from` folders exist
  stopifnot(all(dir.exists(folders$from)))
  
  # Iteratively copy folders & contents
  success <- 
    cl_lapply(split(folders, collapse::seq_row(folders)), 
              .cl = cl, 
              .fun = function(folder) {
                # Clean directory, if required
                # * This is important as old files may inappropriately persist otherwise 
                if (dir.exists(folder$to)) {
                  unlink(folder$to, recursive = TRUE)
                }
                # Create directory, if required
                if (!dir.exists(folder$to)) {
                  dir.create(folder$to, recursive = TRUE)
                } 
                # Copy folder & contents
                file.copy(folder$from, dirname(folder$to), recursive = TRUE, overwrite = TRUE)
              })
  
  # Return success
  folders[, success := unlist(success)]
  folders
  
}

# Example
# from <- c("/Users/lavended/Desktop/1", 
#           "/Users/lavended/Desktop/2")
# to   <- c("/Users/lavended/Desktop/new/1", 
#           "/Users/lavended/Desktop/new/2")
# dirs.copy(from, to, cl = 2L)

#' log.txt

sink_open <- function(log.folder = NULL, log.txt = NULL) {
  # Define log.txt path
  if (!is.null(log.folder)) {
    # Define name
    if (is.null(log.txt)) {
      log.txt <- paste0("log-", as.numeric(Sys.time()), ".txt")
    }
    # Define full file path & validate 
    log.txt <- file.path(log.folder, log.txt)
    # stopifnot(!file.exists(log.txt))
    # Define connection 
    log.txt <- file(log.txt, open = "wt")
    # Open sink
    sink(log.txt, append = TRUE)
    sink(log.txt, type = "message", append = TRUE)
    # Print start time
    print(Sys.time())
  }
  invisible(log.txt)
}

sink_close <- function(log.txt = NULL) {
  if (!is.null(log.txt)) {
    print(Sys.time())
    sink()
    sink(type = "message")
  }
  invisible(NULL)
}

#' Take a coffee break (pause and rest the computer)
#' Only one interval & duration is currently supported 

coffee <- function(computer = "MCC02XT0AZJGH5",
                   folder = here_data("coffee"), 
                   interval = 2 * 60 * 60,
                   duration = 15 * 60
                   ) {
  
  # Take breaks on selected computer(s) only
  if (!(Sys.info()["nodename"] %in% computer)) {
    return(nothing())
  }
  
  # Identify time since last 'coffee' break
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
  time.qs <- file.path(folder, "time.qs")
  if (!file.exists(time.qs)) {
    qs::qsave(Sys.time(), time.qs)
    return(nothing())
  }
  time              <- qs::qread(time.qs)
  time_since_coffee <- difftime(Sys.time(), time, units = "secs")
  
  # Take a break if needed
  if (time_since_coffee >= interval) {
    # Sleep for `duration`
    Sys.sleep(duration)
    # Record the time of the break
    # (We'll sleep in another `duration` secs)
    qs::qsave(Sys.time(), time.qs)
  }
  
  nothing()
  
}

#' Sort bathymetry datasets

sort_bathysets <- function(x) {
  x |> dirname() |> basename() |> unique() |> gtools::mixedsort()
}

fill <- function(x, name) {
  c(name, rep("", length(x) - 1))
}

#' System utilities

on_server <- function() {
  Sys.info()[["nodename"]] == "siam-linux20"
}
