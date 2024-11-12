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
  
  # Handle duplicates
  folders <- 
    data.table(from = from, to = to) |> 
    group_by(from) |> 
    slice(1L) |> 
    as.data.table()
  from <- folders$from
  to   <- folders$to
  
  # Define folder list
  # * Each element is a unique from -> to pairing 
  folders <- mapply(function(f, t) list(from = f, to = t), from, to, SIMPLIFY = FALSE)
  names(folders) <- NULL
  
  # Validate from folders exist
  stopifnot(all(dir.exists(from)))
  
  # Iteratively copy folders & contents
  success <- cl_lapply(folders, .cl = cl, .fun = function(folder) {
    # Create directory, if required
    if (!dir.exists(folder$to)) {
      dir.create(folder$to, recursive = TRUE)
    }
    # Copy folder & contents
    file.copy(folder$from, dirname(file$to), recursive = TRUE, overwrite = TRUE)
  })
  
  # Return success
  data.table(from = from, to = to, success = success)
  
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
    stopifnot(!file.exists(log.txt))
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

#' Sort bathymetry datasets

sort_bathysets <- function(x) {
  x |> dirname() |> basename() |> unique() |> gtools::mixedsort()
}

fill <- function(x, name) {
  c(name, rep("", length(x) - 1))
}
