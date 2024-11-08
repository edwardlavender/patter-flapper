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

#' log.txt

sink_open <- function(log.folder = NULL) {
  log.txt <- NULL
  if (!is.null(log.folder)) {
    # Define file name
    log.txt <- paste0("log-", as.numeric(Sys.time()), ".txt")
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
